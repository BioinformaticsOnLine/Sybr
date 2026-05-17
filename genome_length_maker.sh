#!/usr/bin/env bash
# genome_length_maker.sh
#
# Generate all_sequence_lengths.txt for syntenyPlotteR / Sybr pipeline.
#
# Usage:
#   ./genome_length_maker.sh -r <reference_species> [options]
#
# Required:
#   -r, --reference   Basename of the reference FASTA (without extension)
#
# Optional:
#   -i, --input       Folder containing all FASTA files (default: current directory)
#   -o, --output      Output file path (default: <input_folder>/all_sequence_lengths.txt)
#   -k, --resolutions Comma- or space-separated resolutions in kb (default: 100,300,500)
#   -h, --help        Show this help message
#
# Examples:
#   ./genome_length_maker.sh -r Genus_sps1
#   ./genome_length_maker.sh -r Genus_sps1 -i inputs/fasta
#   ./genome_length_maker.sh -r Genus_sps1 -i inputs/fasta -o inputs/synteny_processing/all_sequence_lengths.txt
#   ./genome_length_maker.sh -r Genus_sps1 -i inputs/fasta -k 100,300,500,1000
#   ./genome_length_maker.sh -r Genus_sps1 -i inputs/fasta -o my_lengths.txt -k 200 400

set -euo pipefail

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

usage() {
    cat <<EOF
Usage: $(basename "$0") -r <reference_species> [options]

Required:
  -r, --reference   Basename of the reference FASTA (no extension)

Optional:
  -i, --input       Folder containing all FASTA files (default: current directory)
  -o, --output      Output file path (default: <input_folder>/all_sequence_lengths.txt)
  -k, --resolutions Comma- or space-separated resolutions in kb (default: 100,300,500)
  -h, --help        Show this help

Examples:
  $(basename "$0") -r Genus_sps1
  $(basename "$0") -r Genus_sps1 -i inputs/fasta
  $(basename "$0") -r Genus_sps1 -i inputs/fasta -o inputs/synteny_processing/all_sequence_lengths.txt
  $(basename "$0") -r Genus_sps1 -i inputs/fasta -k 100,300,500,1000
  $(basename "$0") -r Genus_sps1 -i inputs/fasta -o my_lengths.txt -k 200 400
EOF
    exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────
REFERENCE=""
INPUT_DIR="."
OUTPUT_FILE=""
RESOLUTION_STR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reference)
            REFERENCE="$2"; shift 2 ;;
        -i|--input)
            INPUT_DIR="$2"; shift 2 ;;
        -o|--output)
            OUTPUT_FILE="$2"; shift 2 ;;
        -k|--resolutions)
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do
                RESOLUTION_STR="${RESOLUTION_STR} $1"
                shift
            done ;;
        -h|--help)
            usage ;;
        *)
            error "Unknown argument: $1. Use -h for help." ;;
    esac
done

# ── Validate required args ────────────────────────────────────────────────────
[ -z "$REFERENCE" ] && error "Reference species name is required. Use -r <name>."

# ── Resolve input directory ───────────────────────────────────────────────────
INPUT_DIR="${INPUT_DIR%/}"
[ -d "$INPUT_DIR" ] || error "Input directory not found: $INPUT_DIR"

# ── Resolve output file ───────────────────────────────────────────────────────
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="${INPUT_DIR}/all_sequence_lengths.txt"
fi
mkdir -p "$(dirname "$OUTPUT_FILE")"

# ── Parse resolutions ─────────────────────────────────────────────────────────
if [ -z "$RESOLUTION_STR" ]; then
    RESOLUTIONS=(100 300 500)
else
    IFS=', ' read -r -a raw_res <<< "${RESOLUTION_STR//,/ }"
    RESOLUTIONS=()
    for r in "${raw_res[@]}"; do
        [ -z "$r" ] && continue
        RESOLUTIONS+=("${r%k}")
    done
fi

# ── Validate seqkit ───────────────────────────────────────────────────────────
if ! command -v seqkit &>/dev/null; then
    error "seqkit not found. Install with: conda install -c bioconda seqkit"
fi

# ── Find reference FASTA ──────────────────────────────────────────────────────
REF_FA=""
for ext in fa fasta fna; do
    candidate="${INPUT_DIR}/${REFERENCE}.${ext}"
    if [ -f "$candidate" ]; then
        REF_FA="$candidate"
        break
    fi
done
[ -z "$REF_FA" ] && error "Reference FASTA not found in '${INPUT_DIR}': tried ${REFERENCE}.fa / .fasta / .fna"

info "Input dir  : $INPUT_DIR"
info "Reference  : $REF_FA"
info "Resolutions: ${RESOLUTIONS[*]/%/k}"
info "Output     : $OUTPUT_FILE"

# ── Helper: sorted chr lengths from a FASTA ──────────────────────────────────
fasta_lengths() {
    local fa="$1"
    seqkit fx2tab --name --length "$fa" \
        | awk 'BEGIN{OFS="\t"} {print $1, $2}' \
        | sort -V -k1,1
}

# ── Write output ──────────────────────────────────────────────────────────────
> "$OUTPUT_FILE"

# Reference: one block per resolution
info "Processing reference: $REFERENCE"
ref_lengths=$(fasta_lengths "$REF_FA")

for res in "${RESOLUTIONS[@]}"; do
    label="${REFERENCE}:${res}k"
    while IFS=$'\t' read -r chr len; do
        printf '%s\t%s\t%s\n' "$chr" "$len" "$label"
    done <<< "$ref_lengths"
done >> "$OUTPUT_FILE"

# Query species: every other FASTA in INPUT_DIR
info "Scanning for query genomes in: $INPUT_DIR"
query_count=0

for fa in "${INPUT_DIR}"/*.fa "${INPUT_DIR}"/*.fasta "${INPUT_DIR}"/*.fna; do
    [ -e "$fa" ] || continue
    species="$(basename "${fa%.*}")"
    [ "$species" = "$REFERENCE" ] && continue

    info "  Query: $species  ($fa)"
    while IFS=$'\t' read -r chr len; do
        printf '%s\t%s\t%s\n' "$chr" "$len" "$species"
    done < <(fasta_lengths "$fa") >> "$OUTPUT_FILE"
    query_count=$((query_count + 1))
done

# ── Summary ───────────────────────────────────────────────────────────────────
total_lines=$(wc -l < "$OUTPUT_FILE")
echo ""
info "Output written to : $OUTPUT_FILE"
info "Reference         : $REFERENCE  (${#RESOLUTIONS[@]} resolution blocks)"
info "Query genomes     : $query_count"
info "Total lines       : $total_lines"
