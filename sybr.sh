#!/bin/bash
# sybr1.sh — Snakemake runner with spinner, flags, logging, and clean Ctrl+C handling

set -euo pipefail

#######################################
# Trap for Ctrl+C / termination
#######################################
cleanup() {
    echo >&2
    echo -e "${YELLOW}[INFO] Interrupted. Cleaning up...${NC}" >&2

    # Kill spinner first
    [[ -n "${SPINNER_PID:-}" ]] && kill "$SPINNER_PID" 2>/dev/null || true
    
    # Kill snakemake and all its children
    if [[ -n "${SNAKEMAKE_PID:-}" ]]; then
        # Kill the entire process group
        # Get the process group ID
        PGID=$(ps -o pgid= -p "$SNAKEMAKE_PID" 2>/dev/null | grep -o '[0-9]*' || true)
        
        if [[ -n "$PGID" ]]; then
            # Kill the entire process group (negative PGID kills the whole group)
            kill -TERM -"$PGID" 2>/dev/null || true
            sleep 1
            # Force kill if still running
            kill -KILL -"$PGID" 2>/dev/null || true
        else
            # Fallback: kill the process and its children
            kill -TERM "$SNAKEMAKE_PID" 2>/dev/null || true
            sleep 1
            kill -KILL "$SNAKEMAKE_PID" 2>/dev/null || true
        fi
        
        # Additional cleanup: find and kill any remaining snakemake-related processes
        pkill -f "snakemake" 2>/dev/null || true
        
        # Kill perl processes from this pipeline
        pkill -f "synteny_assign_v3_80_genes.pl" 2>/dev/null || true
        
        # Kill any satsuma processes
        pkill -f "Satsuma" 2>/dev/null || true
        
        # Kill any python processes running pipeline scripts
        pkill -f "validate_satsuma_files.py" 2>/dev/null || true
        
        # Kill any bash processes running our pipeline commands
        pkill -f "/usr/bin/bash -c set -euo pipefail" 2>/dev/null || true
        
        # Kill any processes in the current session that might be orphaned
        pkill -s 0 -f "snakemake" 2>/dev/null || true
    fi

    # Remove temporary files
    rm -f /tmp/snakemake_pid_$$ 2>/dev/null || true
    rm -f /tmp/sybr_snakemake_*.log 2>/dev/null || true
    rm -f /tmp/sybr_merged_*.yaml 2>/dev/null || true
    rm -f synteny_params.yaml 2>/dev/null || true
    
    # Restore cursor
    tput cnorm >&2
    
    # Exit with interrupt code
    exit 130
}
trap cleanup INT TERM

#######################################
# Cleanup on normal exit
#######################################
normal_exit_cleanup() {
    # Remove temporary files on normal exit
    rm -f /tmp/snakemake_pid_$$ 2>/dev/null || true
    rm -f /tmp/sybr_snakemake_*.log 2>/dev/null || true
    rm -f /tmp/sybr_merged_*.yaml 2>/dev/null || true
    rm -f synteny_params.yaml 2>/dev/null || true
}
trap normal_exit_cleanup EXIT

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
get_cpu_count() {
    if command -v nproc > /dev/null; then
        nproc
    elif command -v sysctl > /dev/null; then
        sysctl -n hw.ncpu
    else
        echo 1
    fi
}

CONFIG_FILE="run_sybr_config.yaml"
PATHS_FILE="pipeline_paths.yaml"
CORES=$(get_cpu_count)
TARGET="all"
LOG_FILE=""
UNLOCK=false
DRY_RUN=false
KEEP_GOING=false
VERBOSE=false
SKIP_VALIDATION=false

# New parameters for synteny assignment
WINDOW_SIZES=""  # Empty string means use defaults
STEP_SIZE=""     # Empty string means use defaults
CUSTOM_RESOLUTIONS=false  # Flag to track if custom resolutions were provided

#######################################
# Logging helpers
#######################################
log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

#######################################
# Header / Usage
#######################################
print_header() {
    echo -e "${YELLOW}"
    cat <<'EOF'

┏━┓╻ ╻┏┓ ┏━┓
┗━┓┗┳┛┣┻┓┣┳┛
┗━┛ ╹ ┗━┛╹┗╸
============
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
  -P, --paths FILE       Static paths file (default: pipeline_paths.yaml)
  -j, --cores N          Number of cores (default: all available)
  -t, --target RULE      Target rule (default: all)
  -l, --log FILE         Log output to file
  -u, --unlock           Unlock working directory
  -n, --dry-run          Dry run (simulate pipeline)
  -k, --keep-going       Keep going on independent job failures
  -v, --verbose          Verbose Snakemake output
  -s, --skip-validation  Skip input validation
  -w, --window-sizes     Comma-separated window sizes in bp (e.g., 100000,300000,500000)
  -p, --step-size        Step size in bp for synteny assignment (default: 30000)
  -h, --help             Show this help

Examples:
  $0 -c config.yaml -j 8                         # Run with default settings
  $0 --window-sizes 200000,400000 --step-size 50000  # Custom window sizes and step size
  $0 --window-sizes 100000                        # Single window size
  $0 --step-size 25000                           # Custom step size with default windows

Note: Window sizes and step size only affect synteny_assign rules. If not specified,
      defaults are: window sizes = 100000,300000,500000 and step size = 30000.

EOF
    exit 0
}

#######################################
# Argument parsing
#######################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config) CONFIG_FILE="$2"; shift 2 ;;
        -P|--paths) PATHS_FILE="$2"; shift 2 ;;
        -j|--cores) CORES="$2"; shift 2 ;;
        -t|--target) TARGET="$2"; shift 2 ;;
        -l|--log) LOG_FILE="$2"; shift 2 ;;
        -u|--unlock) UNLOCK=true; shift ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        -k|--keep-going) KEEP_GOING=true; shift ;;
        -s|--skip-validation) SKIP_VALIDATION=true; shift ;;
        -v|--verbose) VERBOSE=true; shift ;;
        -w|--window-sizes) 
            WINDOW_SIZES="$2"
            CUSTOM_RESOLUTIONS=true
            shift 2 
            ;;
        -p|--step-size) 
            STEP_SIZE="$2"
            CUSTOM_RESOLUTIONS=true
            shift 2 
            ;;
        -h|--help) usage ;;
        --) shift; break ;;
        -*) log_error "Invalid argument: $1"; usage ;;
        *) 
            # Capture positional arguments as targets if needed, 
            # though the script seems to prefer -t.
            # For now, just break or handle as error.
            log_error "Unsupported positional argument: $1"; usage
            ;;
    esac
done

#######################################
# Validation
#######################################
[ -f "$CONFIG_FILE" ] || { log_error "Config file not found: $CONFIG_FILE"; exit 1; }
[ -f "$PATHS_FILE"  ] || { log_error "Paths file not found: $PATHS_FILE"; exit 1; }
[[ "$CORES" =~ ^[0-9]+$ && "$CORES" -gt 0 ]] || { log_error "Invalid core count"; exit 1; }

# Validate window sizes if provided
if [ -n "$WINDOW_SIZES" ]; then
    # Remove any whitespace
    WINDOW_SIZES=$(echo "$WINDOW_SIZES" | tr -d '[:space:]')
    
    # Split by comma and validate each is a positive integer
    IFS=',' read -ra WINDOWS <<< "$WINDOW_SIZES"
    for window in "${WINDOWS[@]}"; do
        if ! [[ "$window" =~ ^[0-9]+$ && "$window" -gt 0 ]]; then
            log_error "Invalid window size: $window. Must be a positive integer"
            exit 1
        fi
    done
    
    log_info "Using custom window sizes: $WINDOW_SIZES"
else
    # Use default window sizes
    WINDOW_SIZES="100000,300000,500000"
fi

# Validate step size if provided
if [ -n "$STEP_SIZE" ]; then
    if ! [[ "$STEP_SIZE" =~ ^[0-9]+$ && "$STEP_SIZE" -gt 0 ]]; then
        log_error "Invalid step size: $STEP_SIZE. Must be a positive integer"
        exit 1
    fi
    log_info "Using custom step size: $STEP_SIZE"
else
    # Use default step size
    STEP_SIZE="30000"
    log_info "Using default step size: $STEP_SIZE"
fi

#######################################
# Logging redirection
#######################################
if [ -n "$LOG_FILE" ]; then
    exec > >(tee -a "$LOG_FILE") 2>&1
fi

#######################################
# Get script base directory from config
#######################################
get_script_base() {
    # Merge both config files (paths file first, user config overrides)
    if command -v python3 &> /dev/null; then
        local script_base
        script_base=$(python3 - "$PATHS_FILE" "$CONFIG_FILE" <<'PYEOF2'
import yaml, sys, os

def load(path):
    try:
        with open(path, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}

config = load(sys.argv[1])   # pipeline_paths.yaml
config.update(load(sys.argv[2]))  # run_sybr_config.yaml  (wins on conflict)

scripts = config.get("scripts", "script_base")
if not os.path.isabs(scripts):
    scripts = os.path.join(os.getcwd(), scripts)
print(scripts)
PYEOF2
)
        echo "$script_base"
    else
        echo "script_base"
    fi
}

#######################################
# Create parameter file for Snakemake
#######################################
create_param_file() {
    local param_file="synteny_params.yaml"
    
    cat > "$param_file" <<EOF
# Synteny assignment parameters
# This file is auto-generated by the runner script
# Custom parameters: $([ "$CUSTOM_RESOLUTIONS" = true ] && echo "YES" || echo "NO")

synteny_params:
  # Window sizes in base pairs (comma-separated)
  window_sizes: [$WINDOW_SIZES]
  
  # Step size in base pairs
  step_size: $STEP_SIZE
  
  # Flag to indicate custom parameters were provided
  custom_resolutions: $CUSTOM_RESOLUTIONS
EOF
    
    echo "$param_file"
}

#######################################
# Merge both YAML configs into one temp file
#######################################
merge_configs() {
    python3 - "$PATHS_FILE" "$CONFIG_FILE" <<'PYEOF2'
import yaml, sys, os, tempfile

def load(path):
    try:
        with open(path, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}

def deep_merge(base, override):
    """Merge override into base; sub-dicts are merged key-by-key."""
    result = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = {**result[k], **v}
        else:
            result[k] = v
    return result

def rebase_paths(cfg, in_base, out_base):
    """Replace inputs/ and outputs/ prefixes with the configured base dirs."""
    if isinstance(cfg, dict):
        return {k: rebase_paths(v, in_base, out_base) for k, v in cfg.items()}
    if isinstance(cfg, list):
        return [rebase_paths(v, in_base, out_base) for v in cfg]
    if isinstance(cfg, str) and not os.path.isabs(cfg):
        if cfg.startswith("inputs/") or cfg == "inputs":
            return os.path.join(in_base, cfg[len("inputs/"):])
        if cfg.startswith("outputs/") or cfg == "outputs":
            return os.path.join(out_base, cfg[len("outputs/"):])
    return cfg

# pipeline_paths.yaml is the base; run_sybr_config.yaml wins on conflict
config = deep_merge(load(sys.argv[1]), load(sys.argv[2]))

# Apply base_input_dir / base_output_dir if set
in_base  = config.get("base_input_dir",  "inputs")
out_base = config.get("base_output_dir", "outputs")
if in_base != "inputs" or out_base != "outputs":
    config = rebase_paths(config, in_base, out_base)

tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml",
                                   prefix="sybr_merged_", delete=False)
yaml.dump(config, tmp)
tmp.close()
print(tmp.name)
PYEOF2
}

#######################################
# Input Validation Function
#######################################
validate_inputs() {
    log_info "Validating input files..."

    # Check if Python is available
    if ! command -v python3 &> /dev/null; then
        log_warning "python3 not found, skipping validation"
        return 0
    fi

    # Merge both configs into a single temp file that validate script can read
    MERGED_CONFIG=$(merge_configs)
    if [ -z "$MERGED_CONFIG" ] || [ ! -f "$MERGED_CONFIG" ]; then
        log_warning "Could not merge config files, skipping validation"
        return 0
    fi

    # Get script_base directory from merged configs
    SCRIPT_BASE=$(get_script_base)

    # Check if script_base directory exists
    if [ ! -d "$SCRIPT_BASE" ]; then
        log_warning "Script base directory not found: $SCRIPT_BASE"
        log_warning "Skipping validation."
        rm -f "$MERGED_CONFIG"
        return 0
    fi

    # Check for validation script in script_base
    VALIDATION_SCRIPT="$SCRIPT_BASE/validate_satsuma_files.py"
    if [ ! -f "$VALIDATION_SCRIPT" ]; then
        log_warning "Validation script not found: $VALIDATION_SCRIPT"
        log_warning "Skipping validation. Please ensure validate_satsuma_files.py is in $SCRIPT_BASE"
        rm -f "$MERGED_CONFIG"
        return 0
    fi

    # Run validation against the single merged config file
    if python3 "$VALIDATION_SCRIPT" --config "$MERGED_CONFIG" --verbose; then
        log_success "Input validation passed"
        rm -f "$MERGED_CONFIG"
        return 0
    else
        log_error "Input validation failed"
        rm -f "$MERGED_CONFIG"
        if [ "$KEEP_GOING" = false ]; then
            log_error "Use --skip-validation to bypass validation or fix input files"
            exit 1
        else
            log_warning "Continuing despite validation errors (--keep-going)"
            return 0
        fi
    fi
}

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
    CMD=(snakemake --configfile "$PATHS_FILE" --configfile "$CONFIG_FILE" -j "$CORES" --rerun-incomplete)

    # Add parameters file if custom resolutions were specified
    if [ "$CUSTOM_RESOLUTIONS" = true ]; then
        PARAM_FILE=$(create_param_file)
        CMD+=(--configfiles "$PARAM_FILE")
        log_info "Using custom synteny parameters from: $PARAM_FILE"
    fi

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
    
    # Input validation
    if ! $SKIP_VALIDATION; then
        validate_inputs
    else
        log_warning "Skipping input validation"
    fi
    
    if $UNLOCK; then
        log_info "Unlocking working directory..."
        snakemake --configfile "$PATHS_FILE" --configfile "$CONFIG_FILE" --unlock || true
        echo
    fi

    build_snakemake_cmd

    start_time=$(date +%s)

    ###################################
    # Run Snakemake
    ###################################
    # Temp file captures all snakemake output so we can show it on failure
    SNAKEMAKE_LOG=$(mktemp /tmp/sybr_snakemake_XXXX.log)
    # Don't set trap EXIT here - we have a global trap

    set +e
    # Run snakemake in a process group so we can kill the whole group later
    setsid "${CMD[@]}" > >(
        tee "$SNAKEMAKE_LOG" | while IFS= read -r line; do
            if [[ "$line" == *Completed* ]]; then
                printf "\n${GREEN}✔ %s${NC}\n" "$line"
            fi
        done
    ) 2>&1 &

    SNAKEMAKE_PID=$!
    echo "$SNAKEMAKE_PID" > /tmp/snakemake_pid_$$

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
        log_error "Pipeline failed after ${duration}s — Snakemake output:"
        echo "────────────────────────────────────────────────" >&2
        cat "$SNAKEMAKE_LOG" >&2
        echo "────────────────────────────────────────────────" >&2
        exit 1
    fi
}

main "$@"
