import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description="Robust consensus genome selection per MAG using ORF hits"
)

parser.add_argument("--profiling", type=str, required=True,
                    help="DIAMOND/MMseqs output file (outfmt 6, TSV)")
parser.add_argument("--metadata", type=str, required=True,
                    help="Metadata table (TSV)")
parser.add_argument("--mag_map", type=str, required=True,
                    help="TSV mapping contig -> MAG")
parser.add_argument("--output", type=str, required=True,
                    help="Output file (TSV)")

args = parser.parse_args()


def select_best_genome_per_mag(profiling_file, metadata_file, mag_map_file, output_file):

    # --- 1. Lire résultats DIAMOND/MMseqs ---
    df = pd.read_csv(
        profiling_file,
        sep="\t",
        header=None,
        names=[
            "qseqid", "sseqid", "evalue", "bits", "pident",
            "alnlen", "qcov", "tcov"
        ]
    )

    # conversion numérique
    numeric_cols = ["evalue", "bits", "pident", "alnlen", "qcov"]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)

    # --- 2. Filtrer hits faibles ---
    df = df[(df["evalue"] < 1e-5) & (df["qcov"] > 0.5)]

    if df.empty:
        print("⚠ Aucun hit après filtrage.")
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        return

    # --- 3. Meilleur hit par ORF ---
    df = df.sort_values(
        by=["qseqid", "evalue", "pident", "bits"],
        ascending=[True, True, False, False]
    )

    df_best_orf = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    # --- 4. Extraire contig depuis ORF ---
    df_best_orf["contig"] = df_best_orf["qseqid"].str.rsplit("_", n=1).str[0]

    # --- 5. Charger mapping contig -> MAG ---
    mag_map = pd.read_csv(mag_map_file, sep="\t")

    df_best_orf = df_best_orf.merge(
        mag_map,
        on="contig",
        how="left"
    )

    df_best_orf = df_best_orf.dropna(subset=["mag"])

    if df_best_orf.empty:
        print("Aucun ORF associé à un MAG après merge.")
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        return

    # --- 6. Calcul scores ---
    df_best_orf["evalue"] = df_best_orf["evalue"].replace(0, np.finfo(float).eps)
    df_best_orf["log_score"] = -np.log10(df_best_orf["evalue"])
    df_best_orf["pident_weighted"] = df_best_orf["pident"] * df_best_orf["alnlen"]

    grouped = (
        df_best_orf.groupby(["mag", "sseqid"])
        .agg(
            n_orf=("qseqid", "count"),
            score_sum=("log_score", "sum"),
            pident_weighted_sum=("pident_weighted", "sum"),
            alnlen_sum=("alnlen", "sum"),
            best_evalue=("evalue", "min")
        )
        .reset_index()
    )

    grouped["pident"] = grouped["pident_weighted_sum"] / grouped["alnlen_sum"]
    grouped = grouped.drop(columns=["pident_weighted_sum", "alnlen_sum"])

    # total ORFs par MAG
    total_orf = (
        df_best_orf.groupby("mag")["qseqid"]
        .count()
        .reset_index()
        .rename(columns={"qseqid": "total_orf"})
    )

    grouped = grouped.merge(total_orf, on="mag")
    grouped["prop_orf"] = grouped["n_orf"] / grouped["total_orf"]

    # score final pondéré (harmonisé : intègre pident, comme dans les scripts IMG)
    grouped["score_final"] = grouped["score_sum"] * grouped["prop_orf"] * (grouped["pident"] / 100)

    # --- 7. Sélection du meilleur génome par MAG (tri + drop_duplicates, harmonisé) ---
    grouped = grouped.sort_values(
        by=["mag", "score_final", "n_orf", "best_evalue"],
        ascending=[True, False, False, True]
    )
    best_genomes = grouped.drop_duplicates(subset=["mag"], keep="first").copy()

    best_genomes = best_genomes.rename(columns={"sseqid": "genome_id"})
    best_genomes["genome_id"] = best_genomes["genome_id"].str.rsplit("_", n=1).str[0]

    # --- 8. Ajouter metadata ---
    metadata_df = pd.read_csv(metadata_file, sep="\t")

    merged = best_genomes.merge(
        metadata_df,
        on="genome_id",
        how="left"
    )

    # --- 9. Sauvegarde ---
    merged.to_csv(output_file, sep="\t", index=False)


# --- Exécution ---
select_best_genome_per_mag(
    args.profiling,
    args.metadata,
    args.mag_map,
    args.output
)