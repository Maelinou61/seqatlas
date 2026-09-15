import argparse
import pandas as pd
import numpy as np
import re
import duckdb

parser = argparse.ArgumentParser(
    description="Robust consensus genome selection per contig using MMseqs/DIAMOND results with IMG metadata"
)

parser.add_argument("--profiling", type=str, required=True,
                    help="DIAMOND/MMseqs output file (outfmt 6, TSV)")
parser.add_argument("--metadata", type=str, required=True,
                    help="Metadata table (TSV or parquet)")
parser.add_argument("--taxonomy", type=str, required=True,
                    help="Taxonomy table (TSV or parquet)")
parser.add_argument("--output", type=str, required=True,
                    help="Output file (TSV)")

args = parser.parse_args()


def extract_uvig_id(sseqid):
    """Extrait IMGVR_UViG_XXXXXXXXXX_XXXXXX (premier segment avant |)."""
    first = sseqid.split("|")[0].strip()
    if re.match(r'IMGVR_UViG_\d+_\d+', first):
        return first
    return None


def extract_img_genome_id(sseqid):
    """Extrait l'IMG Genome ID depuis le sseqid"""
    match = re.search(r'IMGVR_UViG_(\d+)_\d+', sseqid)
    if match:
        return match.group(1)
    match = re.search(r'Ga(\d+)', sseqid)
    if match:
        return match.group(1)
    if sseqid.isdigit():
        return sseqid
    return None


def select_best_genome_per_contig(profiling_file, metadata_file, taxonomy_file, output_file):

    # --- 1. Lire les résultats ---
    df = pd.read_csv(
        profiling_file,
        sep="\t",
        header=None,
        names=[
            "qseqid", "sseqid", "evalue", "bits", "pident",
            "alnlen", "qcov", "tcov"
        ]
    )

    df["evalue"] = pd.to_numeric(df["evalue"])
    df["bits"]   = pd.to_numeric(df["bits"])
    df["pident"] = pd.to_numeric(df["pident"])
    df["qcov"]   = pd.to_numeric(df["qcov"])
    df["alnlen"] = pd.to_numeric(df["alnlen"])

    # --- 2. Extraction des identifiants ---
    df["uvig_id"]       = df["sseqid"].apply(extract_uvig_id)        # pour join taxonomie
    df["img_genome_id"] = df["sseqid"].apply(extract_img_genome_id)  # pour join métadonnées
    df = df.dropna(subset=["uvig_id"])

    if df.empty:
        print("Aucun IMG Genome ID extrait.")
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        return

    # --- 3. Filtrer hits faibles ---
    df = df[(df["evalue"] < 1e-5) & (df["qcov"] > 0.5)]

    if df.empty:
        print("Aucun hit après filtrage.")
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        return

    # --- 4. Garder le meilleur hit par ORF ---
    df = df.sort_values(
        by=["qseqid", "evalue", "pident", "bits"],
        ascending=[True, True, False, False]
    )
    df_best_orf = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    # --- 5. Extraire le contig ---
    df_best_orf["contig"] = df_best_orf["qseqid"].str.rsplit("_", n=1).str[0]

    # --- 6. Calcul score robuste ---
    df_best_orf["evalue"] = df_best_orf["evalue"].replace(0, np.finfo(float).eps)
    df_best_orf["log_score"] = -np.log10(df_best_orf["evalue"])
    df_best_orf["pident_weighted"] = df_best_orf["pident"] * df_best_orf["alnlen"]

    grouped = (
        df_best_orf.groupby(["contig", "img_genome_id"])
        .agg(
            n_orf=("qseqid", "count"),
            score_sum=("log_score", "sum"),
            pident_weighted_sum=("pident_weighted", "sum"),
            alnlen_sum=("alnlen", "sum"),
            best_evalue=("evalue", "min")
        )
        .reset_index()
    )

    # UVIG du meilleur hit (par score) pour [contig, img_genome_id] -> utilisé pour le join taxonomie
    # (au lieu de l'UVIG le plus fréquent, pour rester cohérent avec le critère de
    # significativité statistique utilisé partout ailleurs dans le pipeline)
    best_uvig_map = (
        df_best_orf.sort_values(
            by=["contig", "img_genome_id", "log_score", "evalue"],
            ascending=[True, True, False, True]
        )
        .drop_duplicates(subset=["contig", "img_genome_id"], keep="first")
        [["contig", "img_genome_id", "uvig_id"]]
    )
    grouped = grouped.merge(best_uvig_map, on=["contig", "img_genome_id"], how="left")

    grouped["pident"] = grouped["pident_weighted_sum"] / grouped["alnlen_sum"]
    grouped = grouped.drop(columns=["pident_weighted_sum", "alnlen_sum"])

    # Nombre total d'ORFs par contig
    total_orf = (
        df_best_orf.groupby("contig")["qseqid"]
        .count()
        .reset_index()
        .rename(columns={"qseqid": "total_orf"})
    )

    grouped = grouped.merge(total_orf, on="contig")
    grouped["prop_orf"] = grouped["n_orf"] / grouped["total_orf"]

    # Score final pondéré (harmonisé avec le script MAG : intègre pident)
    grouped["score_final"] = grouped["score_sum"] * grouped["prop_orf"] * (grouped["pident"] / 100)

    # --- 7. Sélection du meilleur génome par contig ---
    grouped = grouped.sort_values(
        by=["contig", "score_final", "n_orf", "best_evalue"],
        ascending=[True, False, False, True]
    )
    genome_consensus = grouped.drop_duplicates(subset=["contig"], keep="first")
    genome_consensus = genome_consensus.rename(columns={"img_genome_id": "genome_id"})

    # --- 8. Ajouter métadonnées + taxonomie ---
    genome_consensus["genome_id"] = genome_consensus["genome_id"].astype(str)
    genome_ids = set(genome_consensus["genome_id"])
    uvig_ids   = set(genome_consensus["uvig_id"].dropna())

    con = duckdb.connect()
    con.execute("CREATE TEMP TABLE _genome_ids AS SELECT unnest(?) AS id", [list(genome_ids)])
    con.execute("CREATE TEMP TABLE _uvig_ids   AS SELECT unnest(?) AS id", [list(uvig_ids)])

    # Taxonomy : join sur uvig
    if taxonomy_file.endswith('.parquet'):
        taxonomy_df = con.execute(f"""
            SELECT * FROM read_parquet('{taxonomy_file}')
            WHERE uvig IN (SELECT id FROM _uvig_ids)
        """).df()
    else:
        taxonomy_df = con.execute(f"""
            SELECT * FROM read_csv('{taxonomy_file}', sep='\t', header=true, all_varchar=true)
            WHERE uvig IN (SELECT id FROM _uvig_ids)
        """).df()

    # Metadata : join sur genome_id (= taxon_oid numérique)
    if metadata_file.endswith('.parquet'):
        metadata_df = con.execute(f"""
            SELECT * REPLACE (CAST("IMG Genome ID" AS VARCHAR) AS "IMG Genome ID")
            FROM read_parquet('{metadata_file}')
            WHERE CAST("IMG Genome ID" AS VARCHAR) IN (SELECT id FROM _genome_ids)
        """).df()
    else:
        metadata_df = con.execute(f"""
            SELECT * FROM read_csv('{metadata_file}', sep='\t', header=true, all_varchar=true)
            WHERE "IMG Genome ID" IN (SELECT id FROM _genome_ids)
        """).df()

    con.close()

    merged_taxonomy = genome_consensus.merge(
        taxonomy_df,
        left_on="uvig_id",
        right_on="uvig",
        how="left"
    )

    merged = merged_taxonomy.merge(
        metadata_df,
        left_on="genome_id",
        right_on="IMG Genome ID",
        how="left"
    )

    # --- 9. Sauvegarde ---
    merged.to_csv(output_file, sep="\t", index=False)


select_best_genome_per_contig(args.profiling, args.metadata, args.taxonomy, args.output)