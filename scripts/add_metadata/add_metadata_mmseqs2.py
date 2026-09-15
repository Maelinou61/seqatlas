import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description="Robust consensus genome selection per contig using MMseqs/DIAMOND results"
)
parser.add_argument("--profiling", type=str, required=True,
                    help="DIAMOND/MMseqs output file (outfmt 6, TSV)")
parser.add_argument("--metadata", type=str, required=True,
                    help="Metadata table (TSV)")
parser.add_argument("--output", type=str, required=True,
                    help="Output file (TSV)")
args = parser.parse_args()


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

    # Conversion numérique (important)
    df["evalue"] = pd.to_numeric(df["evalue"])
    df["bits"] = pd.to_numeric(df["bits"])
    df["pident"] = pd.to_numeric(df["pident"])
    df["qcov"] = pd.to_numeric(df["qcov"])
    df["alnlen"] = pd.to_numeric(df["alnlen"])

    # --- 2. Filtrer hits faibles ---
    df = df[(df["evalue"] < 1e-5) & (df["qcov"] > 0.5)]

    if df.empty:
        print("⚠ Aucun hit après filtrage.")
        pd.DataFrame().to_csv(output_file, sep="\t", index=False)
        return

    # --- 3. Garder le meilleur hit par ORF ---
    df = df.sort_values(
        by=["qseqid", "evalue", "pident", "bits"],
        ascending=[True, True, False, False]
    )
    df_best_orf = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    # --- 4. Extraire le contig ---
    df_best_orf["contig"] = df_best_orf["qseqid"].str.rsplit("_", n=1).str[0]

    # --- 5. Calcul score robuste ---
    # Correction appliquée directement sur df_best_orf (et non sur df, qui n'est
    # plus utilisé ensuite) pour éviter les log10(0) = inf.
    df_best_orf["evalue"] = df_best_orf["evalue"].replace(0, np.finfo(float).eps)
    df_best_orf["log_score"] = -np.log10(df_best_orf["evalue"])

    df_best_orf["pident_weighted"] = df_best_orf["pident"] * df_best_orf["alnlen"]

    grouped = (
        df_best_orf.groupby(["contig", "sseqid"])
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
    # Nettoyage colonnes intermédiaires
    grouped = grouped.drop(columns=["pident_weighted_sum", "alnlen_sum"])

    # Nombre total d'ORFs par contig
    total_orf = df_best_orf.groupby("contig")["qseqid"].count().reset_index()
    total_orf.columns = ["contig", "total_orf"]

    grouped = grouped.merge(total_orf, on="contig")
    grouped["prop_orf"] = grouped["n_orf"] / grouped["total_orf"]

    # Score final pondéré (harmonisé : intègre pident)
    grouped["score_final"] = grouped["score_sum"] * grouped["prop_orf"] * (grouped["pident"] / 100)

    # --- 6. Sélection du meilleur génome par contig ---
    # Tri + drop_duplicates (au lieu d'idxmax seul) pour départager les égalités
    # strictes sur score_final via n_orf puis best_evalue.
    grouped = grouped.sort_values(
        by=["contig", "score_final", "n_orf", "best_evalue"],
        ascending=[True, False, False, True]
    )
    best_genomes = grouped.drop_duplicates(subset=["contig"], keep="first")
    best_genomes = best_genomes.rename(columns={"sseqid": "genome_id"})
    best_genomes["genome_id"] = best_genomes["genome_id"].str.rsplit("_", n=1).str[0]

    # --- 7. Ajouter métadonnées ---
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    merged = best_genomes.merge(
        metadata_df,
        left_on="genome_id",
        right_on="genome_id",
        how="left"
    )

    # --- 8. Sauvegarde ---
    merged.to_csv(output_file, sep="\t", index=False)


# --- Exécution ---
select_best_genome_per_contig(args.profiling, args.metadata, args.output)
