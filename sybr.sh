#!/bin/bash
# sybr2.sh — Snakemake runner with spinner, flags, logging, and clean Ctrl+C handling

set -euo pipefail

#######################################
# Trap for Ctrl+C / termination
#######################################
cleanup() {
    echo >&2
    echo -e "${YELLOW}[INFO] Interrupted. Cleaning up...${NC}" >&2

    [[ -n "${SPINNER_PID:-}" ]] && kill "$SPINNER_PID" 2>/dev/null || true
    [[ -n "${SNAKEMAKE_PID:-}" ]] && kill "$SNAKEMAKE_PID" 2>/dev/null || true

    tput cnorm >&2
    exit 130
}
trap cleanup INT TERM

#######################################
# Color codes (non-bold)
#######################################
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
WHITE='\033[0;37m'
Bright_violet='\033[0;95m'
LIGHT_VIOLET='\033[0;95m'
NC='\033[0m'

#######################################
# Defaults
#######################################
CONFIG_FILE="run_sybr_config.yaml"
CORES=$(nproc)
TARGET="all"
LOG_FILE=""
UNLOCK=false
DRY_RUN=false
KEEP_GOING=false
VERBOSE=false

#######################################
# Logging helpers
#######################################
log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1" >&2; }

#######################################
# Header / Usage
#######################################
print_header() {
    echo -e "${YELLOW}"
    cat <<'EOF'
 ░▒▓███████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░  
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
 ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓███████▓▒░░▒▓███████▓▒░  
       ░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
       ░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓███████▓▒░   ░▒▓█▓▒░   ░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
EOF
    echo -e "${NC}"
    echo "========================================"
    echo "    Snakemake Pipeline"
    echo "========================================"
}


usage() {
    print_header
    cat <<EOF

Usage: $0 [OPTIONS]

Options:
  -c, --config FILE      Configuration file (default: run_sybr_config.yaml)
  -j, --cores N          Number of cores (default: all available)
  -t, --target RULE      Target rule (default: all)
  -l, --log FILE         Log output to file
  -u, --unlock           Unlock working directory
  -n, --dry-run          Dry run (simulate pipeline)
  -k, --keep-going       Keep going on independent job failures
  -v, --verbose          Verbose Snakemake output
  -h, --help             Show this help

EOF
    exit 0
}

#######################################
# Argument parsing (your requested flags)
#######################################
PARSED_ARGS=$(getopt -o c:j:t:l:unkvh \
  --long config:,cores:,target:,log:,unlock,dry-run,keep-going,verbose,help \
  -n "$0" -- "$@")

eval set -- "$PARSED_ARGS"

while true; do
    case "$1" in
        -c|--config) CONFIG_FILE="$2"; shift 2 ;;
        -j|--cores) CORES="$2"; shift 2 ;;
        -t|--target) TARGET="$2"; shift 2 ;;
        -l|--log) LOG_FILE="$2"; shift 2 ;;
        -u|--unlock) UNLOCK=true; shift ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        -k|--keep-going) KEEP_GOING=true; shift ;;
        -v|--verbose) VERBOSE=true; shift ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) log_error "Invalid argument"; usage ;;
    esac
done

#######################################
# Validation
#######################################
[ -f "$CONFIG_FILE" ] || { log_error "Config file not found: $CONFIG_FILE"; exit 1; }
[[ "$CORES" =~ ^[0-9]+$ && "$CORES" -gt 0 ]] || { log_error "Invalid core count"; exit 1; }

#######################################
# Logging redirection
#######################################
if [ -n "$LOG_FILE" ]; then
    exec > >(tee -a "$LOG_FILE") 2>&1
fi

#######################################
# Spinner (stderr, PID-controlled)
#######################################
start_spinner() {
    [[ -t 2 ]] || return

    SPINNER_FRAMES=("◐" "◓" "◑" "◒")
    SPINNER_INDEX=0

    tput civis >&2

    (
        while kill -0 "$SNAKEMAKE_PID" 2>/dev/null; do
            frame=${SPINNER_FRAMES[$SPINNER_INDEX]}
            SPINNER_INDEX=$(( (SPINNER_INDEX + 1) % ${#SPINNER_FRAMES[@]} ))
            printf "\r${LIGHT_VIOLET}Sybr is Running... ${YELLOW}%s${NC}" "$frame" >&2
            sleep 0.1
        done
    ) &
    SPINNER_PID=$!
}

stop_spinner() {
    [[ -n "${SPINNER_PID:-}" ]] && wait "$SPINNER_PID" 2>/dev/null || true
    printf "\r\033[K" >&2
    tput cnorm >&2
}

#######################################
# Build Snakemake command (array-safe)
#######################################
build_snakemake_cmd() {
    CMD=(snakemake --configfile "$CONFIG_FILE" -j "$CORES")

    $DRY_RUN && CMD+=(-n)
    $KEEP_GOING && CMD+=(-k)
    ! $VERBOSE && CMD+=(--quiet)

    CMD+=("$TARGET")
}

#######################################
# Main
#######################################
main() {
    print_header
    echo
    
    # REMOVED: Configuration display section
    # log_info "Configuration:"
    # echo "  Config file:  $CONFIG_FILE"
    # echo "  Cores:        $CORES"
    # echo "  Target:       $TARGET"
    # echo "  Unlock:       $UNLOCK"
    # echo "  Dry run:      $DRY_RUN"
    # echo "  Keep going:   $KEEP_GOING"
    # echo "  Verbose:      $VERBOSE"
    # [ -n "$LOG_FILE" ] && echo "  Log file:     $LOG_FILE"
    # echo

    if $UNLOCK; then
        log_info "Unlocking working directory..."
        snakemake --configfile "$CONFIG_FILE" --unlock || true
        echo
    fi

    build_snakemake_cmd

    # REMOVED: "Starting pipeline execution..." message
    # log_info "Starting pipeline execution..."
    # echo
    # echo "Command: ${CMD[*]}"
    # echo

    start_time=$(date +%s)

    ###################################
    # Run Snakemake
    ###################################
    set +e
    "${CMD[@]}" > >(
        while IFS= read -r line; do
            if [[ "$line" == *Completed* ]]; then
                printf "\n${GREEN}✔ %s${NC}\n" "$line"
            fi
        done
    ) 2>&1 &
    SNAKEMAKE_PID=$!

    start_spinner
    wait "$SNAKEMAKE_PID"
    exit_status=$?
    set -e

    stop_spinner

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    echo
    if [ "$exit_status" -eq 0 ]; then
        log_success "Pipeline completed successfully in ${duration}s"
    else
        log_error "Pipeline failed after ${duration}s"
        exit 1
    fi
}

main "$@"
