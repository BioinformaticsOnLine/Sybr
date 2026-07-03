"""
SYBR API — Server configuration.

All settings are loaded from environment variables with sensible defaults.
"""

import os
from pathlib import Path

# ── Base directories ─────────────────────────────────────────────────────────

# Root of the SYBR pipeline (where Snakefile, sybr.sh, etc. live)
SYBR_PIPELINE_DIR = Path(
    os.getenv("SYBR_PIPELINE_DIR", Path(__file__).resolve().parent.parent)
)

# Directory where per-job input/output directories are created
SYBR_JOBS_DIR = Path(
    os.getenv("SYBR_JOBS_DIR", SYBR_PIPELINE_DIR / "jobs")
)

# SQLite database file
SYBR_DB_PATH = Path(
    os.getenv("SYBR_DB_PATH", SYBR_PIPELINE_DIR / "api" / "sybr_api.db")
)

# ── Server settings ─────────────────────────────────────────────────────────

API_HOST = os.getenv("SYBR_API_HOST", "0.0.0.0")
API_PORT = int(os.getenv("SYBR_API_PORT", "8000"))
API_PREFIX = "/api/v1"

# ── Pipeline settings ───────────────────────────────────────────────────────

# Path to sybr.sh runner script
SYBR_RUNNER = SYBR_PIPELINE_DIR / "sybr.sh"

# Path to pipeline_paths.yaml (static paths, rarely changes)
PIPELINE_PATHS_YAML = SYBR_PIPELINE_DIR / "pipeline_paths.yaml"

# Maximum concurrent pipeline jobs (thread pool size)
MAX_CONCURRENT_JOBS = int(os.getenv("SYBR_MAX_JOBS", "2"))

# Default number of cores per job
DEFAULT_CORES = int(os.getenv("SYBR_DEFAULT_CORES", "4"))

# Maximum Snakemake cores allowed per submitted job
MAX_CORES_PER_JOB = int(os.getenv("SYBR_MAX_CORES_PER_JOB", "20"))

# ── File upload settings ────────────────────────────────────────────────────

# Maximum upload size per file (2 GB)
MAX_UPLOAD_SIZE_BYTES = int(os.getenv("SYBR_MAX_UPLOAD_MB", "2048")) * 1024 * 1024

# Allowed file extensions by category
ALLOWED_EXTENSIONS = {
    "fasta": {".fa", ".fasta", ".fna"},
    "satsuma_alignments": {".txt"},
    "sequence_lengths": {".txt"},
    "scaffolds": {".txt"},
    "lastz_alignments": {".axt"},
    "species_info": {".txt"},
    "tree": {".txt"},
    "classification": {".eba"},
    "annotation": {".tsv"},
    "kegg": {".txt"},
    "hgt": {".txt"},
}

# Category → subdirectory mapping inside a job's inputs/ folder
CATEGORY_PATHS = {
    "fasta": "fasta",
    "satsuma_alignments": "synteny_processing/Satsuma_alignments",
    "sequence_lengths": "synteny_processing",
    "scaffolds": "synteny_processing",
    "lastz_alignments": "Ancestor_seq_recunstruction/LastZ_alignments",
    "species_info": "Ancestor_seq_recunstruction",
    "tree": "Ancestor_seq_recunstruction",
    "classification": "eba_analysis",
    "annotation": "enrichment_analysis",
    "kegg": "enrichment_analysis",
    "hgt": "HGTs",
}

# ── Logging ──────────────────────────────────────────────────────────────────

LOG_LEVEL = os.getenv("SYBR_LOG_LEVEL", "INFO")
