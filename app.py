import csv
import os
import shlex
import shutil
import subprocess
import sys
import uuid
from itertools import islice
from pathlib import Path
from typing import Dict, List

from flask import Flask, render_template, request, send_from_directory, url_for
from werkzeug.utils import secure_filename

BASE_DIR = Path(__file__).resolve().parent


def load_local_env(env_path: Path) -> None:
    if not env_path.exists():
        return

    for raw_line in env_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'\"")
        os.environ.setdefault(key, value)


load_local_env(BASE_DIR / ".env")

UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", BASE_DIR / "uploads"))
RESULTS_DIR = Path(os.getenv("RESULTS_DIR", BASE_DIR / "results"))
MAPS_DIR = Path(os.getenv("MAPS_DIR", RESULTS_DIR / "maps"))
MMSEQS_DB = os.getenv("MMSEQS_DB", "")
MMSEQS_DB_DIR = Path(os.getenv("MMSEQS_DB_DIR", BASE_DIR / "databases"))
MMSEQS_BIN = os.getenv("MMSEQS_BIN", "mmseqs")
MMSEQS_TMP_DIR = Path(os.getenv("MMSEQS_TMP_DIR", BASE_DIR / "tmp"))
MMSEQS_THREADS = os.getenv("MMSEQS_THREADS", "").strip()
MMSEQS_EXTRA_ARGS = shlex.split(os.getenv("MMSEQS_EXTRA_ARGS", ""))
PRODIGAL_GV_BIN = os.getenv("PRODIGAL_GV_BIN", "prodigal-gv")
PRODIGAL_GV_MODE = os.getenv("PRODIGAL_GV_MODE", "meta")
PRODIGAL_GV_EXTRA_ARGS = shlex.split(os.getenv("PRODIGAL_GV_EXTRA_ARGS", ""))
POSTPROCESS_SCRIPT = Path(
    os.getenv(
        "POSTPROCESS_SCRIPT",
        BASE_DIR / "scripts" / "add_metadata" / "add_metadata_mmseqs2.py",
    )
)
METADATA_FILE = os.getenv("METADATA_FILE", "")

RESULT_COLUMNS = [
    "query",
    "target",
    "evalue",
    "bits",
    "pident",
    "alnlen",
    "qcov",
    "tcov",
]

UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
MMSEQS_TMP_DIR.mkdir(parents=True, exist_ok=True)
MAPS_DIR.mkdir(parents=True, exist_ok=True)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 128 * 1024 * 1024  # 128 MB


def resolve_project_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return BASE_DIR / path


def path_for_form(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(BASE_DIR))
    except ValueError:
        return str(path.resolve())


def get_database_choices() -> List[Dict[str, str]]:
    db_dir = resolve_project_path(MMSEQS_DB_DIR)
    choices_by_value: Dict[str, str] = {}

    if MMSEQS_DB:
        default_path = resolve_project_path(MMSEQS_DB)
        choices_by_value[path_for_form(default_path)] = default_path.name

    if db_dir.exists():
        for dbtype_file in sorted(db_dir.glob("*.dbtype")):
            db_path = dbtype_file.with_suffix("")
            if db_path.name.endswith("_h") or not db_path.exists():
                continue
            choices_by_value[path_for_form(db_path)] = db_path.name

    return [
        {"value": value, "label": label}
        for value, label in sorted(choices_by_value.items(), key=lambda item: item[1].lower())
    ]


def get_default_database_value() -> str:
    choices = get_database_choices()
    if MMSEQS_DB:
        default_value = path_for_form(resolve_project_path(MMSEQS_DB))
        if any(choice["value"] == default_value for choice in choices):
            return default_value
    return choices[0]["value"] if choices else ""


def validate_selected_database(raw_value: str) -> str:
    selected_value = raw_value or get_default_database_value()
    allowed_values = {choice["value"] for choice in get_database_choices()}
    if not selected_value or selected_value not in allowed_values:
        raise ValueError("Base MMseqs invalide ou non disponible")
    return selected_value


@app.context_processor
def inject_database_choices() -> Dict[str, object]:
    selected_db = request.form.get("mmseqs_db") or get_default_database_value()
    return {
        "db_choices": get_database_choices(),
        "selected_db": selected_db,
    }


def allowed_file(filename: str) -> bool:
    lower = filename.lower()
    return lower.endswith((".faa", ".fa", ".fasta", ".fna", ".ffn"))


def looks_like_nucleotide_fasta(path: Path) -> bool:
    if path.suffix.lower() == ".faa":
        return False

    sequence_chars = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_chars.extend(char.upper() for char in line if not char.isspace())

    if not sequence_chars:
        raise ValueError("Aucune sequence exploitable dans le fichier")

    nucleotide_chars = set("ACGTUNRYSWKMBDHV.-")
    invalid_chars = {char for char in sequence_chars if char not in nucleotide_chars}
    return not invalid_chars


def normalize_sequence_input(raw_text: str) -> str:
    text = (raw_text or "").strip()
    if not text:
        raise ValueError("Aucune sequence fournie")

    if text.startswith(">"):
        return text if text.endswith("\n") else f"{text}\n"

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    sequence = "".join(lines).replace(" ", "")
    if not sequence:
        raise ValueError("La sequence est vide")

    allowed_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")
    invalid = sorted({char for char in sequence.upper() if char not in allowed_chars})
    if invalid:
        invalid_preview = "".join(invalid[:10])
        raise ValueError(f"Sequence invalide, caracteres non supportes: {invalid_preview}")

    return f">query_sequence\n{sequence}\n"


def run_mmseqs(query_fasta: Path, output_tsv: Path, mmseqs_db: str) -> subprocess.CompletedProcess:
    if not mmseqs_db:
        raise ValueError("MMSEQS_DB n'est pas configuré")

    mmseqs_path = shutil.which(MMSEQS_BIN) if not os.path.sep in MMSEQS_BIN else MMSEQS_BIN
    if not mmseqs_path or not Path(mmseqs_path).exists():
        raise FileNotFoundError(
            f"Binaire MMseqs introuvable: '{MMSEQS_BIN}'. "
            "Installe mmseqs2 ou définis MMSEQS_BIN avec le chemin complet."
        )

    cmd = [
        mmseqs_path,
        "easy-search",
        str(query_fasta),
        str(resolve_project_path(mmseqs_db)),
        str(output_tsv),
        str(MMSEQS_TMP_DIR),
        "-s", "7.5", 
        "-c", "0.5", 
        "--max-seqs", "50",
        "--format-output",
        ",".join(RESULT_COLUMNS),
    ]

    if MMSEQS_THREADS:
        cmd.extend(["--threads", MMSEQS_THREADS])

    return subprocess.run(cmd, capture_output=True, text=True, check=False)   


def run_prodigal_gv(input_fasta: Path, output_faa: Path) -> subprocess.CompletedProcess:
    prodigal_path = (
        shutil.which(PRODIGAL_GV_BIN)
        if os.path.sep not in PRODIGAL_GV_BIN
        else PRODIGAL_GV_BIN
    )
    if not prodigal_path or not Path(prodigal_path).exists():
        raise FileNotFoundError(
            f"Binaire prodigal-gv introuvable: '{PRODIGAL_GV_BIN}'. "
            "Installe prodigal-gv ou définis PRODIGAL_GV_BIN avec le chemin complet."
        )

    output_gff = output_faa.with_suffix(".gff")
    cmd = [
        prodigal_path,
        "-i",
        str(input_fasta),
        "-a",
        str(output_faa),
        "-o",
        str(output_gff),
        "-p",
        PRODIGAL_GV_MODE,
    ]
    cmd.extend(PRODIGAL_GV_EXTRA_ARGS)

    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def read_result_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []

    rows: List[Dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, fieldnames=RESULT_COLUMNS, delimiter="\t")
        for row in reader:
            rows.append(row)

    return rows


def read_headered_tsv_rows(path: Path, limit: int | None = None) -> List[Dict[str, str]]:
    if not path.exists():
        return []

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if limit is None:
            return list(reader)
        return list(islice(reader, limit))


def run_metadata_postprocess(input_tsv: Path, output_tsv: Path) -> subprocess.CompletedProcess:
    if not POSTPROCESS_SCRIPT.exists():
        raise FileNotFoundError(f"Script de post-traitement introuvable: {POSTPROCESS_SCRIPT}")
    if not METADATA_FILE:
        raise ValueError("METADATA_FILE n'est pas configuré")

    cmd = [
        sys.executable,
        str(POSTPROCESS_SCRIPT),
        "--profiling",
        str(input_tsv),
        "--metadata",
        METADATA_FILE,
        "--output",
        str(output_tsv),
    ]

    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def generate_kepler_map(input_tsv: Path, output_html: Path) -> int:
    import pandas as pd
    from keplergl import KeplerGl

    df = pd.read_csv(input_tsv, sep="\t")
    if df.empty or "latitude" not in df.columns or "longitude" not in df.columns:
        return 0

    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    points = df.dropna(subset=["latitude", "longitude"]).copy()
    if points.empty:
        return 0

    group_cols = ["latitude", "longitude"]
    if "biosample_name" in points.columns:
        group_cols.insert(0, "biosample_name")

    agg_spec = {
        "n_hits": ("genome_id", "size"),
    }
    if "genome_id" in points.columns:
        agg_spec["n_genomes"] = ("genome_id", "nunique")
    if "study_name" in points.columns:
        agg_spec["study_name"] = ("study_name", "first")
    if "microntology" in points.columns:
        agg_spec["microntology"] = ("microntology", "first")
    if "geographic_location" in points.columns:
        agg_spec["geographic_location"] = ("geographic_location", "first")

    points = (
        points.groupby(group_cols, dropna=False)
        .agg(**agg_spec)
        .reset_index()
    )

    if "biosample_name" not in points.columns:
        points["biosample_name"] = "unknown"
    if "n_genomes" not in points.columns:
        points["n_genomes"] = points["n_hits"]

    points = points.sort_values(
        by=["n_hits", "n_genomes"],
        ascending=[False, False],
    )

    if points.empty:
        return 0

    if len(points) > 3000:
        points = points.head(3000).copy()

    config = {
        "version": "v1",
        "config": {
            "visState": {
                "layers": [
                    {
                        "type": "point",
                        "config": {
                            "dataId": "matched_genomes",
                            "label": "Matched genomes",
                            "columns": {
                                "lat": "latitude",
                                "lng": "longitude",
                            },
                            "visConfig": {
                                "opacity": 0.9,
                                "radius": 22,
                                "filled": True,
                                "strokeColor": [255, 255, 255],
                                "strokeOpacity": 0.8,
                                "thickness": 2,
                                "color": [13, 148, 136],
                            },
                        },
                        "visualChannels": {
                            "sizeField": {"name": "n_hits", "type": "integer"},
                            "sizeScale": "sqrt",
                        },
                    }
                ],
                "interactionConfig": {
                    "tooltip": {
                        "enabled": True,
                        "fieldsToShow": {
                            "matched_genomes": [
                                {"name": "biosample_name"},
                                {"name": "n_hits"},
                                {"name": "n_genomes"},
                                {"name": "study_name"},
                                {"name": "microntology"},
                                {"name": "geographic_location"},
                                {"name": "latitude"},
                                {"name": "longitude"},
                            ]
                        },
                    }
                },
            },
            "mapState": {
                "bearing": 0,
                "dragRotate": False,
                "latitude": 20,
                "longitude": 0,
                "pitch": 0,
                "zoom": 1.2,
            },
            "mapStyle": {
                "styleType": "dark-matter-nolabels",
                "visibleLayerGroups": {
                    "label": True,
                    "road": False,
                    "border": True,
                    "building": False,
                    "water": True,
                    "land": True,
                    "3d building": False,
                },
            },
        },
    }

    map_ = KeplerGl(height=700, data={"matched_genomes": points}, config=config)
    map_.save_to_html(file_name=str(output_html), read_only=False)
    return len(points)


@app.get("/")
def index():
    return render_template(
        "index.html",
        rows=[],
        headers=[],
        message="",
        error="",
        map_url="",
        map_points=0,
    )


@app.get("/results/<path:filename>")
def download_result(filename: str):
    return send_from_directory(RESULTS_DIR, filename)


@app.post("/search")
def search():
    uploaded = request.files.get("faa_file")
    sequence_text = (request.form.get("sequence_text") or "").strip()
    selected_db = request.form.get("mmseqs_db", "")
    has_upload = bool(uploaded and uploaded.filename)
    has_sequence = bool(sequence_text)

    try:
        selected_db = validate_selected_database(selected_db)
    except Exception as exc:
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error=str(exc),
            map_url="",
            map_points=0,
        )

    if not has_upload and not has_sequence:
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error="Ajoute un fichier FASTA ou colle une sequence dans la zone de texte",
            map_url="",
            map_points=0,
        )

    if has_upload and not allowed_file(uploaded.filename):
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error="Format non supporté (utilise .faa, .fa, .fasta, .fna ou .ffn)",
            map_url="",
            map_points=0,
        )

    run_id = uuid.uuid4().hex[:12]
    if has_sequence:
        query_file = UPLOAD_DIR / f"{run_id}_pasted_sequence.fasta"
        try:
            query_file.write_text(normalize_sequence_input(sequence_text), encoding="utf-8")
        except Exception as exc:
            return render_template(
                "index.html",
                rows=[],
                headers=[],
                message="",
                error=f"Sequence invalide: {exc}",
                map_url="",
                map_points=0,
            )
    else:
        original_name = secure_filename(uploaded.filename)
        query_file = UPLOAD_DIR / f"{run_id}_{original_name}"
        uploaded.save(query_file)

    result_file = RESULTS_DIR / f"{run_id}_mmseqs.tsv"
    merged_file = RESULTS_DIR / f"{run_id}_metadata.tsv"
    map_file = MAPS_DIR / f"{run_id}_kepler.html"
    mmseqs_query_file = query_file
    message_parts = []
    warnings = []

    try:
        if looks_like_nucleotide_fasta(query_file):
            translated_file = UPLOAD_DIR / f"{run_id}_prodigal_gv.faa"
            completed = run_prodigal_gv(query_file, translated_file)
            if completed.returncode != 0:
                stderr = (completed.stderr or completed.stdout or "").strip()
                raise RuntimeError(stderr or "prodigal-gv a échoué")
            if not translated_file.exists() or translated_file.stat().st_size == 0:
                raise RuntimeError("prodigal-gv n'a produit aucune sequence proteique")

            mmseqs_query_file = translated_file
            message_parts.append(
                f"Sequence nucleotidique detectee, ORFs predites avec prodigal-gv: {translated_file.name}."
            )
    except Exception as exc:
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error=f"Erreur prodigal-gv: {exc}",
            map_url="",
            map_points=0,
        )

    try:
        completed = run_mmseqs(mmseqs_query_file, result_file, selected_db)
    except Exception as exc:
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error=f"Erreur backend: {exc}",
            map_url="",
            map_points=0,
        )

    if completed.returncode != 0:
        stderr = (completed.stderr or "").strip()
        return render_template(
            "index.html",
            rows=[],
            headers=[],
            message="",
            error=f"mmseqs a échoué: {stderr or 'Erreur inconnue'}",
            map_url="",
            map_points=0,
        )

    rows: List[Dict[str, str]]
    map_url = ""
    map_points = 0

    try:
        completed = run_metadata_postprocess(result_file, merged_file)
        if completed.returncode != 0:
            stderr = (completed.stderr or completed.stdout or "").strip()
            raise RuntimeError(stderr or "Le script de post-traitement a echoue")

        rows = read_headered_tsv_rows(merged_file, limit=200)
        try:
            map_points = generate_kepler_map(merged_file, map_file)
            if map_points:
                map_url = url_for("download_result", filename=f"maps/{map_file.name}")
            else:
                warnings.append(
                    "Aucune coordonnee exploitable n'a ete trouvee pour afficher la carte."
                )
        except Exception as exc:
            warnings.append(f"Carte Kepler ignoree: {exc}")
    except Exception as exc:
        rows = read_result_rows(result_file)[:200]
        warnings.append(f"Post-traitement metadata ignore: {exc}")

    headers = list(rows[0].keys()) if rows else RESULT_COLUMNS
    message_parts.append(f"Recherche terminee. Apercu de {len(rows)} ligne(s).")
    message_parts.append(f"Base MMseqs utilisee: {Path(selected_db).name}.")
    message_parts.append(f"Resultat MMseqs brut: {result_file.name}.")
    if merged_file.exists():
        message_parts.append(f"Resultat enrichi: {merged_file.name}.")
    if map_points:
        message_parts.append(f"Carte Kepler generee avec {map_points} point(s).")
    if warnings:
        message_parts.extend(warnings)

    return render_template(
        "index.html",
        rows=rows,
        headers=headers,
        message=" ".join(message_parts),
        error="",
        map_url=map_url,
        map_points=map_points,
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.getenv("PORT", "8000")), debug=False)
