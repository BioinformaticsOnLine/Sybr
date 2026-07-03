# Configuration handling module

import os
import glob
import sys
import subprocess
import tempfile

# Load static paths first, then user config.
# Snakemake does a SHALLOW merge of top-level keys, so dict-valued keys like
# "eba:" that appear in both files would be overwritten entirely.
# We fix this below with an explicit deep-merge after both files are loaded.
_PATHS_CONFIG_FILE = os.getenv("SYBR_PIPELINE_PATHS_CONFIG", "pipeline_paths.yaml")
_RUN_CONFIG_FILE = os.getenv("SYBR_RUN_CONFIG", "run_sybr_config.yaml")

configfile: _PATHS_CONFIG_FILE
configfile: _RUN_CONFIG_FILE

# Deep-merge any top-level dicts that are intentionally split across both files.
# pipeline_paths.yaml is loaded first so run_sybr_config.yaml wins on conflict.
def _deep_merge(base, override):
    """Recursively merge override into base dict, returning a new dict.
    Nested dicts are merged recursively; all other types are overwritten.
    """
    result = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = v
    return result

def _deep_merge_from_files(cfg, *yaml_files):
    """Re-read each yaml file and deep-merge its dict keys into cfg."""
    import yaml as _yaml
    layers = []
    for path in yaml_files:
        try:
            with open(path) as _f:
                layers.append(_yaml.safe_load(_f) or {})
        except Exception:
            layers.append({})
    for layer in layers:
        cfg = _deep_merge(cfg, layer)
    return cfg

config = _deep_merge_from_files(config, _PATHS_CONFIG_FILE, _RUN_CONFIG_FILE)

# Rebase all "inputs/..." and "outputs/..." paths in the config using
# base_input_dir and base_output_dir from run_sybr_config.yaml.
# This lets users set a single pair of absolute dirs and have every
# downstream path automatically point to the right place.
def _rebase_paths(cfg, in_base, out_base):
    """
    Walk every string value in cfg (recursively).
    - Values starting with "inputs/"  → replace that prefix with in_base
    - Values starting with "outputs/" → replace that prefix with out_base
    Already-absolute paths are left untouched.
    """
    if isinstance(cfg, dict):
        return {k: _rebase_paths(v, in_base, out_base) for k, v in cfg.items()}
    if isinstance(cfg, list):
        return [_rebase_paths(v, in_base, out_base) for v in cfg]
    if isinstance(cfg, str) and not os.path.isabs(cfg):
        if cfg.startswith("inputs/") or cfg == "inputs":
            return os.path.join(in_base, cfg[len("inputs/"):])
        if cfg.startswith("outputs/") or cfg == "outputs":
            return os.path.join(out_base, cfg[len("outputs/"):])
    return cfg

_in_base  = config.get("base_input_dir",  "inputs")
_out_base = config.get("base_output_dir", "outputs")

# Only rebase if explicit base dirs were provided (i.e. not the bare defaults)
if _in_base != "inputs" or _out_base != "outputs":
    config = _rebase_paths(config, _in_base, _out_base)

# ── Determine whether the de-novo alignment module is active ────────────────
_RUN_SATSUMA_ALIGNMENT = config.get("run_stages", {}).get("run_satsuma_alignment", False)

try:
    # When run_satsuma_alignment is ON, satsuma_alignments points to the
    # *output* of SatsumaSynteny2, so it won't exist yet at config time — that
    # is expected and not an error.  We therefore only require the key to be
    # present (it still needs to be defined as a path for the sort rule), but
    # we do not require the directory to already exist.
    required_keys = [
        "synteny_results", "scripts", "reference_name", "reference_species", "run_stages"
    ]
    # satsuma_alignments is always required as a config key (used in sort rule),
    # but when run_satsuma_alignment is active its value is overridden below.
    if not _RUN_SATSUMA_ALIGNMENT:
        required_keys.append("satsuma_alignments")

    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # When de-novo alignment is active, satsuma_align section must be present
    if _RUN_SATSUMA_ALIGNMENT:
        if "satsuma_align" not in config:
            raise ValueError(
                "run_stages.run_satsuma_alignment is true but 'satsuma_align' "
                "section is missing from pipeline_paths.yaml"
            )
        for _sub in ("fasta_dir", "output_dir", "tool_path"):
            if _sub not in config["satsuma_align"]:
                raise ValueError(
                    f"Missing required key satsuma_align.{_sub} in pipeline_paths.yaml"
                )

    # When de-novo LastZ alignment is active, lastz_align section must be present.
    # Also requires satsuma_align.fasta_dir since FASTA files are shared.
    _RUN_LASTZ_ALIGNMENT = config.get("run_stages", {}).get("run_lastz_alignment", False)
    if _RUN_LASTZ_ALIGNMENT:
        if "lastz_align" not in config:
            raise ValueError(
                "run_stages.run_lastz_alignment is true but 'lastz_align' "
                "section is missing from pipeline_paths.yaml"
            )
        if "tool_path" not in config["lastz_align"]:
            raise ValueError(
                "Missing required key lastz_align.tool_path in pipeline_paths.yaml"
            )
        if "satsuma_align" not in config or "fasta_dir" not in config.get("satsuma_align", {}):
            raise ValueError(
                "run_stages.run_lastz_alignment is true but satsuma_align.fasta_dir "
                "is missing from pipeline_paths.yaml — lastz_alignment.smk reads "
                "FASTA files from the same directory as satsuma_alignment.smk"
            )

    # Handle optional scaffolds_file
    if "eba_format" in config:
        if "scaffolds_file" not in config["eba_format"]:
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"
        elif config["eba_format"]["scaffolds_file"] == "Chromosomes":
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"

except Exception as e:
    raise ValueError(f"Error loading config: {str(e)}")

# Function to resolve relative paths
def resolve_path(path, base_dir=None):
    """Resolve relative paths to absolute paths relative to workflow directory"""
    if path is None:
        return None
    if isinstance(path, str):
        if os.path.isabs(path):
            return path
        if path in ["ALL_CHROMOSOMES", "Chromosomes"]:
            return path
        if base_dir:
            return os.path.join(base_dir, path)
        return os.path.join(workflow.basedir, path)
    return path

def resolve_config_paths(config_dict, base_key=None):
    """Recursively resolve all paths in config dictionary"""
    resolved_config = {}
    for key, value in config_dict.items():
        if isinstance(value, dict):
            resolved_config[key] = resolve_config_paths(value, base_key=key)
        elif isinstance(value, str) and any(path_term in key.lower() for path_term in ['dir', 'file', 'path', 'results', 'alignments', 'input']):
            if base_key == 'eba_tools' and key == 'eba_script_path':
                resolved_config[key] = resolve_path(value, workflow.basedir)
            else:
                resolved_config[key] = resolve_path(value)
        else:
            resolved_config[key] = value
    return resolved_config

# Path resolution is handled by _rebase_paths above using base_input_dir /
# base_output_dir.  The resolve_config_paths pass below is intentionally
# removed because it re-anchored every path to workflow.basedir, which
# silently overrode the user-supplied base_input_dir / base_output_dir values.


# ── Global variables accessible in all modules ──────────────────────────────

# When run_satsuma_alignment is active, SATSUMA_ALIGNMENTS points to the
# output_dir of the alignment module so the sort rule can find the .txt files
# produced by SatsumaSynteny2.  Otherwise it uses the static pre-computed dir.
if _RUN_SATSUMA_ALIGNMENT:
    SATSUMA_ALIGNMENTS = config["satsuma_align"]["output_dir"]
    # Ensure the output directory exists so downstream glob won't fail before
    # the alignment rules have created it.
    try:
        os.makedirs(SATSUMA_ALIGNMENTS, exist_ok=True)
    except PermissionError:
        raise ValueError(
            f"[config] Cannot create Satsuma alignment output directory: {SATSUMA_ALIGNMENTS}\n"
            f"Check that base_output_dir in run_sybr_config.yaml points to a writable location."
        )
    print(
        f"[config] run_satsuma_alignment = ON  "
        f"→ alignments will be written to: {SATSUMA_ALIGNMENTS}"
    )
else:
    SATSUMA_ALIGNMENTS = config["satsuma_alignments"]
    print(
        f"[config] run_satsuma_alignment = OFF "
        f"→ using pre-computed alignments from: {SATSUMA_ALIGNMENTS}"
    )

# ── LastZ alignment dir — redirect when run_lastz_alignment is active ────────
_RUN_LASTZ_ALIGNMENT = config.get("run_stages", {}).get("run_lastz_alignment", False)
if _RUN_LASTZ_ALIGNMENT:
    LASTZ_ALIGNMENTS = config["lastz_align"]["output_dir"]
    try:
        os.makedirs(LASTZ_ALIGNMENTS, exist_ok=True)
    except PermissionError:
        raise ValueError(
            f"[config] Cannot create LastZ alignment output directory: {LASTZ_ALIGNMENTS}\n"
            f"Check that base_output_dir in run_sybr_config.yaml points to a writable location."
        )
    # Redirect chainNet.lastZ_alignments so alignment_processing.smk reads
    # from the correct location automatically
    config["chainNet"]["lastZ_alignments"] = LASTZ_ALIGNMENTS
    print(
        f"[config] run_lastz_alignment  = ON  "
        f"→ .axt files will be written to: {LASTZ_ALIGNMENTS}"
    )
else:
    LASTZ_ALIGNMENTS = config["chainNet"]["lastZ_alignments"]
    print(
        f"[config] run_lastz_alignment  = OFF"
        f" → using pre-computed .axt files from: {LASTZ_ALIGNMENTS}"
    )

SYNTENY_RESULTS = config["synteny_results"]
SCRIPTS         = config["scripts"]
REFERENCE       = config["reference_species"]
RUN_STAGES      = config["run_stages"]

# ── Chr-FASTA preparation ─────────────────────────────────────────────────────
# Modules run_satsuma_alignment, chainNet_generation, and
# Ancestor_seq_recunstruction need FASTA files with:
#   • filename  : genus prefix stripped  → sps1.fa  (not Genus_sps1.fa)
#   • headers   : plain numbers prefixed → >chr1    (not >1)
# The original fasta_dir is kept intact; converted copies go to chr_fasta_dir.
# chr-FASTA (short filenames + >chr* headers) is only needed by the
# chainNet / LastZ / Ancestor modules.  Satsuma, synteny_processing,
# eba_analysis, and enrichment_analysis use the original FASTA files
# (Genus_sps1.fa with plain numeric headers like >1, >2 …).
_NEEDS_CHR_FASTA = (
    RUN_STAGES.get("run_lastz_alignment", False)
    or RUN_STAGES.get("chainNet_generation", False)
    or RUN_STAGES.get("Ancestor_seq_recunstruction", False)
)

# chr_fasta_dir: read from config (set in pipeline_paths.yaml).
# The value in pipeline_paths.yaml uses "outputs/fasta_chr" so _rebase_paths
# will have already mapped it to base_output_dir/fasta_chr.
# Fallback (safety net): place it inside base_output_dir, never inside inputs/.
_raw_chr_fasta_dir = config.get("satsuma_align", {}).get(
    "chr_fasta_dir",
    os.path.join(_out_base, "fasta_chr")
)
CHR_FASTA_DIR = _raw_chr_fasta_dir

if _NEEDS_CHR_FASTA:
    try:
        os.makedirs(CHR_FASTA_DIR, exist_ok=True)
    except PermissionError:
        raise ValueError(
            f"[config] Cannot create chr-FASTA directory: {CHR_FASTA_DIR}\n"
            f"Check that base_input_dir in run_sybr_config.yaml points to a writable location."
        )
    print(
        f"[config] chr-FASTA preparation = ON\n"
        f"  source : {config.get('satsuma_align', {}).get('fasta_dir', 'inputs/fasta')}\n"
        f"  dest   : {CHR_FASTA_DIR}"
    )

    # ── Rule: convert one FASTA file ─────────────────────────────────────────
    # Input:  Genus_sps1.fa  with headers like  >1  >2  …
    # Output: sps1.fa        with headers like  >chr1  >chr2  …
    #
    # Wildcard {sps_stem} = species-only part of the filename (e.g. "sps1").
    # The rule discovers the source file by globbing for *_{sps_stem}.fa[sta].
    # For the reference (e.g. Genus_sps1.fa) sps_stem == reference_species.

    _ORIG_FASTA_DIR = config["satsuma_align"]["fasta_dir"]

    # Build a map:  sps_stem (short) → original FASTA path
    # Convention: the short name is everything after the first underscore,
    # lowercased.  If no underscore, the full stem is used as-is.
    #   Genus_sps1.fa  → sps1
    #   sps1.fa        → sps1   (already short, nothing to strip)
    _orig_fasta_map = {}   # sps_stem → original path
    for _ext in (".fa", ".fna", ".fasta"):
        for _f in glob.glob(os.path.join(_ORIG_FASTA_DIR, f"*{_ext}")):
            _full_stem = os.path.splitext(os.path.basename(_f))[0]
            _short     = _full_stem.split("_", 1)[-1]   # "Genus_sps1" → "sps1"
            _orig_fasta_map[_short] = _f

    _CHR_FASTA_SPECIES = sorted(_orig_fasta_map.keys())

    rule prepare_chr_fasta:
        """Copy one FASTA from fasta_dir to chr_fasta_dir.

        Transformations applied:
          1. Filename: strip Genus_ prefix  →  sps1.fa
          2. Headers:  prepend 'chr' to bare numeric IDs  →  >chr1, >chr2 …
             Headers that already start with 'chr' are left unchanged.
        """
        input:
            orig = lambda wildcards: _orig_fasta_map[wildcards.sps_stem]
        output:
            chr_fa = CHR_FASTA_DIR + "/{sps_stem}.fa"
        log:
            CHR_FASTA_DIR + "/logs/{sps_stem}_prepare.log"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {log})"
            echo "=== prepare_chr_fasta ===" | tee {log}
            echo "Source : {input.orig}"     | tee -a {log}
            echo "Dest   : {output.chr_fa}"  | tee -a {log}

            # Rewrite headers: ">1" → ">chr1", ">chr1" left as-is
            awk '
                /^>/ {{
                    header = substr($0, 2)          # strip leading >
                    if (header !~ /^chr/) {{
                        # bare numeric or any non-chr name: prepend chr
                        print ">chr" header
                    }} else {{
                        print $0
                    }}
                    next
                }}
                {{ print }}
            ' "{input.orig}" > "{output.chr_fa}"

            NSEQ=$(grep -c "^>" "{output.chr_fa}" || true)
            echo "Headers rewritten ($NSEQ sequences)." | tee -a {log}
            echo -e "\e[32;1m✔\e[0m \e[33;1mPrepared\e[0m chr-FASTA: \e[36;1m{wildcards.sps_stem}\e[0m" >&2
            """

    # ── Aggregate marker: all chr-FASTAs ready ───────────────────────────────
    rule prepare_chr_fasta_all:
        input:
            expand(CHR_FASTA_DIR + "/{sps_stem}.fa", sps_stem=_CHR_FASTA_SPECIES)
        output:
            touch(CHR_FASTA_DIR + "/.chr_fasta_done")
        shell:
            """
            echo "All chr-FASTA files prepared."
            ls -lh {input}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m all chr-FASTA preparations" >&2
            """

    # ── Cleanup: remove fasta_chr once all consumers have finished ────────────
    # The .fa files in fasta_chr/ are read only by fasta_to_2bit and
    # generate_chrom_sizes, which produce .2bit and .size files respectively.
    # Everything downstream (axt_to_chain, chainPreNet, chainNet, netSyntenic,
    # reorganize_data, deschrambler) uses only those derived files, so it is
    # safe to delete fasta_chr/ as soon as all .2bit and .size files exist.
    # Only defined when chainNet_generation is active (it owns those rules).
    if RUN_STAGES.get("chainNet_generation", False):
        _CHAIN_OUT = config["chainNet"].get("output_dir", ".")

        rule cleanup_chr_fasta:
            """Delete fasta_chr/ after all consumers have finished.

            Waits for:
              • all .2bit files  (fasta_to_2bit)
              • all .size files  (generate_chrom_sizes)
              • lastz_alignment_all marker (run_lastz reads .fa files directly)
              • satsuma_alignment_all marker (run_satsuma reads .fa files directly)
            Only the markers that correspond to active modules are included.
            """
            input:
                twobit = expand(_CHAIN_OUT + "/2bit/{species}.2bit",
                                species=_CHR_FASTA_SPECIES),
                sizes  = expand(_CHAIN_OUT + "/fasize/{species}.size",
                                species=_CHR_FASTA_SPECIES),
                # Wait for LastZ / Satsuma alignments to finish — those rules
                # read the .fa files directly and must complete before deletion.
                # run_lastz reads .fa files directly; wait for it to finish
                lastz_done = (
                    [os.path.join(config["lastz_align"]["output_dir"], ".lastz_alignment_done")]
                    if RUN_STAGES.get("run_lastz_alignment", False)
                    else []
                )
                # Note: run_satsuma_alignment uses the ORIGINAL fasta_dir, not
                # CHR_FASTA_DIR, so no dependency on satsuma_done needed here.
            output:
                touch(_CHAIN_OUT + "/.chr_fasta_cleanup_done")
            params:
                chr_fasta_dir = CHR_FASTA_DIR
            shell:
                r"""
                set -euo pipefail
                echo -e "\e[33;1mCleaning up\e[0m chr-FASTA directory: {params.chr_fasta_dir}" >&2

                # Safety: confirm all derived files exist before deleting anything
                MISSING=0
                for f in {input.twobit} {input.sizes}; do
                    if [ ! -f "$f" ]; then
                        echo "  ERROR: expected derived file missing: $f" >&2
                        MISSING=1
                    fi
                done
                if [ "$MISSING" -eq 1 ]; then
                    echo "  Aborting cleanup — some derived files are missing." >&2
                    exit 1
                fi

                rm -rf "{params.chr_fasta_dir}"
                echo -e "\e[32;1m✔\e[0m \e[33;1mRemoved\e[0m chr-FASTA directory: {params.chr_fasta_dir}" >&2
                """

# Define ebrs_dir globally
ebrs_dir = config["eba_format"].get("ebrs_dir", "EBRs")

# Get synteny parameters if provided
if "synteny_params" in config:
    WINDOW_SIZES       = config["synteny_params"]["window_sizes"]
    STEP_SIZE          = config["synteny_params"]["step_size"]
    CUSTOM_RESOLUTIONS = config["synteny_params"].get("custom_resolutions", False)
else:
    WINDOW_SIZES       = [100000, 300000, 500000]
    STEP_SIZE          = 30000
    CUSTOM_RESOLUTIONS = False

# Create output directory if it doesn't exist
try:
    os.makedirs(SYNTENY_RESULTS, exist_ok=True)
except PermissionError:
    raise ValueError(
        f"[config] Cannot create synteny results directory: {SYNTENY_RESULTS}\n"
        f"Check that base_output_dir in run_sybr_config.yaml points to a writable location.\n"
        f"Current value: {config.get('base_output_dir', '(not set)')}"
    )

# ── Pre-flight input validation ──────────────────────────────────────────────
# Calls validate_satsuma_files.py with the fully merged + rebased config so
# every required input file/directory is checked before the DAG is built.
def _run_input_validation(cfg):
    validator = os.path.join(workflow.basedir, "script_base", "validate_satsuma_files.py")
    if not os.path.isfile(validator):
        print(f"[config] WARNING: Input validator not found at {validator} — skipping pre-flight check")
        return

    try:
        import yaml as _yaml_val
    except ImportError:
        print("[config] WARNING: PyYAML not available — skipping pre-flight input check")
        return

    # Dump the merged + rebased config to a temp file so the validator can read it
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as _tmp:
        _yaml_val.dump(dict(cfg), _tmp)
        _tmp_path = _tmp.name

    try:
        result = subprocess.run(
            [sys.executable, validator, "--config", _tmp_path]
        )
        if result.returncode != 0:
            raise ValueError(
                "[config] Pre-flight input validation failed.\n"
                "Resolve the missing files / directories listed above, then re-run."
            )
    finally:
        try:
            os.unlink(_tmp_path)
        except OSError:
            pass

_run_input_validation(config)

# ── Sample discovery ─────────────────────────────────────────────────────────
if RUN_STAGES.get("synteny_processing", True):
    try:
        if _RUN_SATSUMA_ALIGNMENT:
            # Samples are derived from the FASTA files in fasta_dir.
            # The alignment module names its output files as
            # {query}_{reference_name}.txt, so we discover from FASTA stems.
            _fasta_dir = config["satsuma_align"]["fasta_dir"]
            _ref_name  = config["reference_name"]
            _fasta_stems = []
            for _ext in (".fa", ".fna", ".fasta"):
                _fasta_stems.extend(
                    os.path.splitext(os.path.basename(f))[0]
                    for f in glob.glob(os.path.join(_fasta_dir, f"*{_ext}"))
                )
            # Query species = all stems except the reference
            _query_stems = sorted(s for s in _fasta_stems if s != _ref_name)
            # Sample name = just the query species name, matching the
            # {query}.txt output file from run_satsuma
            SAMPLES = _query_stems
            if not SAMPLES:
                raise ValueError(
                    f"[config] run_satsuma_alignment is ON but no query FASTA "
                    f"files found in {_fasta_dir} (excluding reference '{_ref_name}')"
                )
            print(f"[config] Samples derived from FASTA files: {SAMPLES}")
        else:
            # Original behaviour: discover from existing alignment .txt files
            sample_files = glob.glob(f"{SATSUMA_ALIGNMENTS}/*_*.txt")
            SAMPLES = sorted(
                os.path.splitext(os.path.basename(f))[0]
                for f in sample_files
            )
            if not SAMPLES:
                raise ValueError(
                    f"No sample files found in {SATSUMA_ALIGNMENTS}"
                )
            print(f"Discovered {len(SAMPLES)} samples: {SAMPLES}")

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error discovering samples: {str(e)}")
else:
    SAMPLES = []

# ── Conditional inputs for the top-level 'all' rule ─────────────────────────
def get_conditional_inputs():
    inputs = []

    # NOTE: CHR_FASTA_DIR/.chr_fasta_done is intentionally NOT listed here.
    # It is an intermediate marker used only for rule-to-rule dependency
    # chaining (via chr_ready inputs in fasta_to_2bit / run_satsuma /
    # run_lastz).  The cleanup rule deletes the whole fasta_chr/ directory
    # (including that marker), so listing it here would cause a MissingFiles
    # error in the final 'all' rule.  The cleanup marker below already proves
    # the full preparation + consumption chain completed successfully.

    # If de-novo alignment is active, require its completion marker first
    if RUN_STAGES.get("run_satsuma_alignment", False):
        _align_out = config["satsuma_align"]["output_dir"]
        inputs.append(os.path.join(_align_out, ".satsuma_alignment_done"))

    # Synteny processing outputs (only if enabled)
    if RUN_STAGES.get("synteny_processing", True) and SAMPLES:
        inputs.extend([
            expand(
                f"{SYNTENY_RESULTS}/{{sample}}_out/{{window_size}}/st_input",
                sample=SAMPLES,
                window_size=WINDOW_SIZES
            ),
        ])

        # Add outputs for each window size dynamically
        for window_size in WINDOW_SIZES:
            inputs.extend(
                expand(
                    f"{SYNTENY_RESULTS}/synteny_results/{window_size}/{{sample}}_out/synteny_assign_done",
                    sample=SAMPLES
                )
            )

        # Add EBA format and overlap resolution output markers
        if config.get("out_final"):
            inputs.extend([
                os.path.join(SYNTENY_RESULTS, config["out_final"], ".done_marker"),
                os.path.join(SYNTENY_RESULTS, config["out_final"], ".overlap_resolve_done"),
                os.path.join(SYNTENY_RESULTS, ".concatination_done"),
                os.path.join(SYNTENY_RESULTS, ".reformatting_done"),
                os.path.join(SYNTENY_RESULTS, ".plots_generated_done"),
                os.path.join(SYNTENY_RESULTS, ".split_by_species_done"),
                os.path.join(SYNTENY_RESULTS, ".linear_synteny_plots_done")
            ])

    # EBA analysis outputs (only if enabled)
    if RUN_STAGES.get("eba_analysis", True):
        if config.get("eba_format", {}).get("pre_EBA_dir"):
            pre_EBA_dir = os.path.join(workflow.basedir, config["eba_format"]["pre_EBA_dir"])
            inputs.extend([
                os.path.join(pre_EBA_dir, "EBA-input", ".copy_complete"),
                os.path.join(pre_EBA_dir, "EBA_WORKING/eba_complete.marker"),
                os.path.join(pre_EBA_dir, "EBA_WORKING/EBA_OUT/Merge/Result_Merge.final")
            ])

        if config.get("eba_format", {}).get("mshsbs_dir"):
            mshsbs_dir = os.path.join(workflow.basedir, config["eba_format"]["mshsbs_dir"])
            inputs.extend([
                os.path.join(mshsbs_dir, "msHSBs.txt"),
                os.path.join(mshsbs_dir, "msHSBs_done.txt"),
                os.path.join(mshsbs_dir, "msHSBs_sequences.fasta"),
                os.path.join(mshsbs_dir, "msHSBs_fasta_done.txt"),
            ])

        if config.get("eba_format", {}).get("ebrs_dir"):
            ebrs_dir = os.path.join(workflow.basedir, config["eba_format"]["ebrs_dir"])
            inputs.extend([
                os.path.join(ebrs_dir, "EBRs.txt"),
                os.path.join(ebrs_dir, "ebrs_processed.marker"),
                os.path.join(ebrs_dir, "EBRs_stats.txt"),
                os.path.join(ebrs_dir, "EBRs_stats.html"),
            ])

    # Enrichment analysis outputs (only if enabled)
    if RUN_STAGES.get("enrichment_analysis", True):
        msHSBs_dir = config["enrichment"].get("msHSBs_dir", config["eba_format"].get("mshsbs_dir", "msHSBs"))
        msHSBs_dir = os.path.join(workflow.basedir, msHSBs_dir)

        enrichment_inputs = [
            os.path.join(msHSBs_dir, "msHSBs_NCBI_genes_overlap.txt"),
            os.path.join(msHSBs_dir, "foreground.txt"),
            os.path.join(msHSBs_dir, "background.txt"),
        ]

##        if config["getenrich"]["r"] == "ko":
##            enrichment_inputs.append(os.path.join(msHSBs_dir, "3kegg_annotationTOgenes.txt"))

        if RUN_STAGES.get("eba_analysis", True) and config.get("eba_format", {}).get("ebrs_dir"):
            ebrs_dir = os.path.join(workflow.basedir, config["eba_format"]["ebrs_dir"])
            all_ebrs_dir = os.path.join(ebrs_dir, "EBRs_Enrichment_results")
            enrichment_inputs.extend([
                os.path.join(all_ebrs_dir, "EBRs_NCBI_genes_overlap.txt"),
                os.path.join(all_ebrs_dir, "foreground.txt"),
                os.path.join(all_ebrs_dir, "background.txt"),
                os.path.join(ebrs_dir, "EBRs_split_done.txt"),
                os.path.join(ebrs_dir, "EBRs_subdirs_overlap_done.txt"),
                
#                expand(
#                    os.path.join(ebrs_dir, "{lineage}_lineage", "getENRICH_done.txt"),
#                    lineage=[
#                        os.path.basename(d).replace("_lineage", "")
#                        for d in glob.glob(f"{ebrs_dir}/*_lineage")
#                    ]
#                ),
            ])

##            if config["getenrich"]["r"] == "ko":
##                enrichment_inputs.append(os.path.join(all_ebrs_dir, "3kegg_annotationTOgenes.txt"))

        inputs.extend(enrichment_inputs)

    # If de-novo LastZ alignment is active, require its completion marker
    # before chainNet_generation rules start.
    if RUN_STAGES.get("run_lastz_alignment", False):
        inputs.append(os.path.join(LASTZ_ALIGNMENTS, ".lastz_alignment_done"))

    # Alignment processing outputs (only if enabled)
    if RUN_STAGES.get("chainNet_generation", True):
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            expand(f"{ALIGN_OUTPUT}/2bit/{{species}}.2bit", species=get_fasta_species()),
            f"{ALIGN_OUTPUT}/clean_axt.done",
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}.chain", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}-chainSplit", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/fasize/{{species}}.size", species=get_fasta_species()),
            f"{ALIGN_OUTPUT}/chainPreNet/processing_done.txt",
            f"{ALIGN_OUTPUT}/netFiles.done",
            f"{ALIGN_OUTPUT}/netSyntenic.done",
            f"{ALIGN_OUTPUT}/data/{REFERENCE}/reorganized.done"
        ])

    # Deschrambler outputs (only if enabled)
    if RUN_STAGES.get("Ancestor_seq_recunstruction", True):
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            f"{ALIGN_OUTPUT}/deschrambler.done"
        ])
    # HGT overlap outputs — each sub-analysis is gated on its upstream stage
    if RUN_STAGES.get("hgt_overlap_analysis", True):
        _hgt_out_dir = config.get("hgt_overlap", {}).get(
            "output_dir",
            os.path.join(config.get("base_output_dir", "outputs"), "hgt_overlap")
        )
        # msHSBs and EBRs overlaps require eba_analysis
        if RUN_STAGES.get("eba_analysis", True):
            inputs.extend([
                os.path.join(_hgt_out_dir, "hgt_overlapping_msHSBs.txt"),
                os.path.join(_hgt_out_dir, "hgt_msHSBs_overlap_summary.txt"),
                os.path.join(_hgt_out_dir, "hgt_overlapping_EBRs.txt"),
                os.path.join(_hgt_out_dir, "hgt_EBRs_overlap_summary.txt"),
            ])
        # ancestoral_genes overlap requires Ancestor_seq_recunstruction (deschrambler)
        if RUN_STAGES.get("Ancestor_seq_recunstruction", True):
            inputs.extend([
                os.path.join(_hgt_out_dir, "hgt_overlapping_ancestoral_genes.txt"),
                os.path.join(_hgt_out_dir, "hgt_ancestoral_genes_overlap_summary.txt"),
            ])
        # Interactive HTML tree requires both eba_analysis and Ancestor_seq_recunstruction
        if RUN_STAGES.get("eba_analysis", True) and RUN_STAGES.get("Ancestor_seq_recunstruction", True):
            inputs.append(os.path.join(_hgt_out_dir, "hgt_tree.html"))


    # chr-FASTA cleanup (only when chainNet_generation is active and chr-FASTA
    # was prepared — cleanup_chr_fasta deletes fasta_chr/ after all .2bit and
    # .size files exist, so it must come after deschrambler in this list to
    # ensure Snakemake sees it as a required final target)
    if _NEEDS_CHR_FASTA and RUN_STAGES.get("chainNet_generation", False):
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.append(f"{ALIGN_OUTPUT}/.chr_fasta_cleanup_done")

    return inputs

# Helper functions for other modules
def get_fasta_species():
    """Get species list for fasta_to_2bit and fasize rules.

    When run_lastz_alignment is ON, .axt files don't exist at DAG-build time
    so globbing LASTZ_ALIGNMENTS returns [].  In that case we derive the query
    species list from the original fasta_dir (always present) using the same
    short-name convention as lastz_alignment.smk, then add the reference.

    When run_lastz_alignment is OFF, pre-computed .axt files exist and we use
    their stems as before.
    """
    if not RUN_STAGES.get("chainNet_generation", True):
        return []

    if RUN_STAGES.get("run_lastz_alignment", False):
        # Derive from original FASTA stems — .axt files don't exist yet
        _orig_dir = config.get("satsuma_align", {}).get("fasta_dir", "")
        _stems = []
        for _ext in (".fa", ".fna", ".fasta"):
            for _f in glob.glob(os.path.join(_orig_dir, f"*{_ext}")):
                _full = os.path.splitext(os.path.basename(_f))[0]
                _short = _full.split("_", 1)[-1]
                _stems.append(_short)
        return sorted(set(_stems))   # includes reference short name
    else:
        axt_files = glob.glob(f"{LASTZ_ALIGNMENTS}/*.axt")
        query_species = [
            os.path.splitext(os.path.basename(f))[0]
            for f in axt_files
            if os.path.splitext(os.path.basename(f))[0] != REFERENCE
        ]
        ref_short = REFERENCE
        return sorted(set(query_species + [ref_short]))

def get_alignment_species():
    """Get list of query species (excludes reference).

    When run_lastz_alignment is ON, derives from original fasta_dir stems
    (same logic as get_fasta_species but excludes the reference).
    When OFF, uses existing .axt file stems.
    """
    if not RUN_STAGES.get("chainNet_generation", True):
        return []

    if RUN_STAGES.get("run_lastz_alignment", False):
        _orig_dir = config.get("satsuma_align", {}).get("fasta_dir", "")
        _stems = []
        for _ext in (".fa", ".fna", ".fasta"):
            for _f in glob.glob(os.path.join(_orig_dir, f"*{_ext}")):
                _full = os.path.splitext(os.path.basename(_f))[0]
                _short = _full.split("_", 1)[-1]
                if _short != REFERENCE:
                    _stems.append(_short)
        return sorted(set(_stems))
    else:
        axt_files = glob.glob(f"{LASTZ_ALIGNMENTS}/*.axt")
        return [
            os.path.splitext(os.path.basename(f))[0]
            for f in axt_files
            if os.path.splitext(os.path.basename(f))[0] != REFERENCE
        ]

# Conditional rule function
def should_run_rule(condition_key):
    """Check if a rule should run based on configuration"""
    return RUN_STAGES.get(condition_key, True)
