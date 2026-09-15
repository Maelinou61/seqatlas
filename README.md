# MMseqs Simple Webserver

Serveur web minimal pour:
- uploader un fichier de séquences (`.faa`, `.fa`, `.fasta`)
- exécuter `mmseqs easy-search` sur une base existante
- enrichir les hits selon la base sélectionnée : VIRE ou MetaVR
- afficher les génomes retenus sur une carte monde Kepler

## 1) Installation

Avec conda, recommandé pour installer aussi `mmseqs2` et `prodigal-gv`:

```bash
cd /home/ubuntu/mmseqs-web
conda env create -f environment.yml
conda activate mmseqs-web
```

Alternative si `mmseqs2` et `prodigal-gv` sont déjà installés séparément:

```bash
cd /home/ubuntu/mmseqs-web
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 2) Configuration

```bash
cp .env.example .env
```

Puis adapte au minimum:
- `MMSEQS_DB` vers ta base mmseqs existante
- `MMSEQS_DB_DIR` vers le dossier contenant les bases proposées dans l'interface
- `MMSEQS_BIN` si `mmseqs` n'est pas dans ton `PATH`
- `PRODIGAL_GV_BIN` si `prodigal-gv` n'est pas dans ton `PATH`
- `MMSEQS_THREADS` pour contrôler le parallélisme
- les variables VIRE ou MetaVR de la section suivante

### Profils de post-traitement

Le serveur choisit le post-traitement à partir de la base MMseqs sélectionnée :

- une base dont le nom contient `metavr` ou `imgvr` utilise le profil **MetaVR** ;
- les autres bases utilisent le profil **VIRE** ;
- si une base MetaVR a un autre nom, ajoute son nom ou son chemin à `METAVR_DATABASES` (valeurs séparées par des virgules).

VIRE conserve la configuration historique :

```dotenv
VIRE_POSTPROCESS_SCRIPT=scripts/add_metadata/add_metadata_mmseqs2.py
VIRE_METADATA_FILE=/chemin/vers/vire_with_coordinates.parquet
```

`POSTPROCESS_SCRIPT` et `METADATA_FILE` restent acceptés pour éviter de casser une configuration existante. `VIRE_*` est prioritaire lorsqu'il est défini.

MetaVR nécessite obligatoirement deux tables distinctes : `--metadata` (métadonnées IMG) et `--taxonomy` (taxonomie UViG). Configure-les ainsi :

```dotenv
METAVR_POSTPROCESS_SCRIPT=scripts/add_metadata/add_metadata_genomes_metavr.py
METAVR_METADATA_FILE=/chemin/vers/Source_dataset_metadata.tsv
METAVR_TAXONOMY_FILE=/chemin/vers/IMGVR5_UViG.tsv

# Facultatif si le nom de la base ne contient pas metavr ou imgvr
METAVR_DATABASES=ma_base_imgvr
```

Le script `add_metadata_genomes_metavr.py` relie la taxonomie au meilleur `uvig_id` et les métadonnées au `IMG Genome ID`. Il accepte les fichiers TSV et Parquet. Les colonnes IMG usuelles (`Latitude`, `Longitude`, `Genome Name / Sample Name`, etc.) sont également adaptées pour la carte Kepler.

Exemple de commande équivalente à ce que lance le serveur:

```bash
mmseqs easy-search \
  uploads/test.faa \
  databases/testDB \
  results/test_mmseqs.tsv \
  tmp \
  --format-output "query,target,evalue,bits,pident,alnlen,qcov,tcov" \
  --threads 8
```

Si `mmseqs` n'est pas dans le `PATH`, définis par exemple:

```bash
MMSEQS_BIN=/chemin/complet/vers/mmseqs
```

## 3) Lancer le serveur

```bash
cd /home/ubuntu/mmseqs-web
conda activate mmseqs-web
set -a; source .env; set +a
python3 app.py
```

Interface: `http://localhost:8000`

## Résultats

- Résultat mmseqs brut: `results/<run_id>_mmseqs.tsv`
- Résultat enrichi metadata: `results/<run_id>_metadata.tsv`
- Carte HTML Kepler: `results/maps/<run_id>_kepler.html`
- L'interface affiche un aperçu des 200 premières lignes.

## Notes

- Le script VIRE `scripts/add_metadata/add_metadata_mmseqs2.py` sélectionne le meilleur génome par contig, puis ajoute les métadonnées `genome_id`.
- Le script MetaVR `scripts/add_metadata/add_metadata_genomes_metavr.py` sélectionne le meilleur génome par contig et ajoute à la fois taxonomie et métadonnées IMG.
- Le fichier `vire_with_coordinates.parquet` contient déjà `genome_id`, `biosample_name`, `latitude`, `longitude` et les métadonnées utiles.
- Si le post-traitement ou la carte échouent, la recherche MMseqs continue quand même et l'interface affiche un avertissement.
- La base `databases/testDB*` déjà présente dans ce dépôt semble être une base MMseqs utilisable telle quelle.
- L'interface liste automatiquement les bases MMseqs trouvées dans `MMSEQS_DB_DIR` via leurs fichiers `.dbtype`.
