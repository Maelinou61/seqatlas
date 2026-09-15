import pandas as pd
from keplergl import KeplerGl

arcs = pd.read_csv("results_bins.tsv", sep="\t")
points = pd.read_csv("biosample_genomes.tsv", sep="\t")

map_ = KeplerGl(height=600)

arcs["identity_scaled"] = (arcs["identity"] - arcs["identity"].min()) / (
    arcs["identity"].max() - arcs["identity"].min()
)

# accentuation visuelle forte
arcs["identity_scaled"] = arcs["identity_scaled"] ** 2

# si n_samples absent → compter automatiquement
if "n_samples" not in points.columns:
    counts = arcs.groupby("biosample_name").size()
    points["n_samples"] = points["biosample_name"].map(counts).fillna(0)

# créer la carte AVEC config
config = {
  "version": "v1",
  "config": {
    "visState": {
      "layers": [

        {
          "type": "arc",
          "config": {
            "dataId": "genome_links",
            "label": "Genome links",
            "columns": {
              "lat0": "latitude_source",
              "lng0": "longitude_source",
              "lat1": "latitude_target",
              "lng1": "longitude_target"
            },
            "visConfig": {
              "opacity": 0.55
            },
            "visualChannels": {
              "sizeField": {"name": "identity_scaled", "type": "real"},
              "sizeScale": "log"
            }
          }
        },

        {
          "type": "point",
          "config": {
            "dataId": "biosamples",
            "label": "Samples",
            "columns": {
              "lat": "latitude_target",
              "lng": "longitude_target"
            }
          },
          "visualChannels": {
            "sizeField": {"name": "n_samples", "type": "integer"},
            "sizeScale": "sqrt"
          }
        }

      ]
    }
  }
}

map_ = KeplerGl(height=700, config=config)

map_.add_data(arcs, "genome_links")
map_.add_data(points, "biosamples")

map_.save_to_html(file_name="genome_contigs_map.html")

print("Carte générée : genome_map.html")
