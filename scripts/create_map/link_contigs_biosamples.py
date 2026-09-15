import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Create a csv file who make the link between the genomes and the biosamples associated")
parser.add_argument("--contig_metadata", type=str, required=True, help="Contig metadata file (tsv)")
parser.add_argument("--geobiosample", type=str, required=True, help="File with the biosample and the latitude - longitude data associated")
parser.add_argument("--output", type=str, required=True, help="Output tsv file with the name of the contig, his coordinates, the target coordinates, the confidence and the environnment")

args = parser.parse_args()

def parse_tables(contig_metadata, geobiosample, output):
    contig_metadata["identity"] = contig_metadata["pident"]
    contig_metadata["longitude_source"] = 45.2889
    contig_metadata["latitude_source"] = -12.7708

    contig_metadata = contig_metadata[["qseqid", "longitude_source", "latitude_source","identity", "biosample_name", "microntology"]]
    contig_metadata = contig_metadata.drop_duplicates(subset=["qseqid"])
    counts = contig_metadata.groupby("biosample_name").size().reset_index(name="n_samples")

    geobiosample = geobiosample[["accession", "longitude", "latitude"]]
    geobiosample = geobiosample.rename(columns={"accession": "biosample_name", "longitude": "longitude_target", "latitude": "latitude_target"})
    
    print(counts)
    merged = contig_metadata.merge(geobiosample, on="biosample_name", how="left")
    merged = merged.merge(counts, on="biosample_name", how="left")
    merged["n_samples"] = merged["n_samples"].fillna(0)

    contigs = merged[["qseqid", "longitude_source", "latitude_source"]]
    biosamples = merged[["biosample_name", "longitude_target", "latitude_target", "microntology","n_samples"]]

    arcs = merged[["qseqid", "longitude_source", "latitude_source", "longitude_target", "latitude_target", "identity"]]

    biosamples = biosamples.copy()
    biosamples.loc[biosamples["microntology"].str.contains('high salinity', na=False), "microntology"] = 'high salinity'
    biosamples.loc[biosamples["microntology"].str.contains('low salinity', na=False), "microntology"] = 'low salinity'

    contigs.to_csv("contigs.tsv", sep="\t", index=False)
    biosamples.to_csv("biosample.tsv", sep="\t", index=False)
    arcs = arcs.drop_duplicates(subset=["qseqid"])
    arcs.to_csv(output, sep="\t", index=False)

def main():
    contig_metadata = pd.read_csv(args.contig_metadata, sep="\t")
    geobiosample = pd.read_csv(args.geobiosample, sep=",")
    output = args.output
    parse_tables(contig_metadata, geobiosample, output)

if __name__ == "__main__":
    main()
