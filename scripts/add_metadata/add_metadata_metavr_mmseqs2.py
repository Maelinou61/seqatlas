import argparse
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser(
    description="Robust consensus genome selection per contig using MMseqs/DIAMOND results with IMG metadata"
)

parser.add_argument("--profiling", type=str, required=True,
                    help="DIAMOND/MMseqs output file (outfmt 6, TSV)")
parser.add_argument("--metadata", type=str, required=True,
                    help="Metadata table (TSV)")
parser.add_argument("--output", type=str, required=True,
                    help="Output file (TSV)")

args = parser.parse_args()


def extract_img_genome_id(sseqid):
    """Extrait l'IMG Genome ID depuis le sseqid"""
    # Pattern principal: IMGVR_UViG_XXXXXXXXX_XXXXXX
    match = re.search(r'IMGVR_UViG_(\d+)_\d+', sseqid)
    if match:
        return match.group(1)
    
    # Alternative: Ga suivi de chiffres
    match = re.search(r'Ga(\d+)', sseqid)
    if match:
        return match.group(1)
    
    # Si c'est déjà juste un ID numérique
    if sseqid.isdigit():
        return sseqid
    
    return None


def select_best_genome_per_contig(profiling_file, metadata_file, output_file):

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

    # Conversion numérique
    df["evalue"] = pd.to_numeric(df["evalue"])
    df["bits"] = pd.to_numeric(df["bits"])
    df["pident"] = pd.to_numeric(df["pident"])
    df["qcov"] = pd.to_numeric(df["qcov"])
    df["alnlen"] = pd.to_numeric(df["alnlen"])

    # --- 2. Extraction de l'IMG Genome ID ---
    df["img_genome_id"] = df["sseqid"].apply(extract_img_genome_id)
    df = df.dropna(subset=["img_genome_id"])

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

    # Calcul du pident pondéré
    grouped["pident"] = grouped["pident_weighted_sum"] / grouped["alnlen_sum"]
    grouped = grouped.drop(columns=["pident_weighted_sum", "alnlen_sum"])

    # Nombre total d'ORFs par contig
    total_orf = df_best_orf.groupby("contig")["qseqid"].count().reset_index()
    total_orf.columns = ["contig", "total_orf"]

    grouped = grouped.merge(total_orf, on="contig")
    grouped["prop_orf"] = grouped["n_orf"] / grouped["total_orf"]

    # Score final pondéré
    grouped["score_final"] = grouped["score_sum"] * grouped["prop_orf"]

    # --- 7. Sélection du meilleur génome par contig ---
    best_idx = grouped.groupby("contig")["score_final"].idxmax()
    best_genomes = grouped.loc[best_idx].rename(columns={"img_genome_id": "genome_id"})

    # --- 8. Ajouter métadonnées ---
    if metadata_file.endswith('.parquet'):
        metadata_df = pd.read_parquet(metadata_file)
    else:
        metadata_df = pd.read_csv(metadata_file, sep="\t")
    
    # Convertir les IMG Genome ID en string pour le merge
    best_genomes["genome_id"] = best_genomes["genome_id"].astype(str)
    metadata_df["IMG Genome ID"] = metadata_df["IMG Genome ID"].astype(str)
    
    merged = best_genomes.merge(
        metadata_df,
        left_on="genome_id",
        right_on="IMG Genome ID",
        how="left"
    )

    # --- 9. Sauvegarde ---
    merged.to_csv(output_file, sep="\t", index=False)


# --- Exécution ---
select_best_genome_per_contig(args.profiling, args.metadata, args.output)
