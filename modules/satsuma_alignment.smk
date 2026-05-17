# ─────────────────────────────────────────────────────────────────────────────
#  satsuma_alignment.smk  —  Optional de-novo Satsuma alignment module
#
#  Activated by setting  run_stages.run_satsuma_alignment: true  in
#  run_sybr_config.yaml.
# ─────────────────────────────────────────────────────────────────────────────

import os
import glob as _glob

if should_run_rule("run_satsuma_alignment"):

    # ── Resolve paths from config ────────────────────────────────────────────
    # Satsuma uses the ORIGINAL fasta_dir with full filenames (Genus_sps1.fa)
    # and original headers (>1, >2 …).  chr-FASTA conversion is NOT needed here.
    _FASTA_DIR  = config["satsuma_align"]["fasta_dir"]
    _OUTPUT_DIR = config["satsuma_align"]["output_dir"]
    _TOOL       = config["satsuma_align"]["tool_path"]
    _THREADS     = config["satsuma_align"].get("threads", 40)
    _MAX_PARALLEL = config["satsuma_align"].get("max_parallel", 1)
    _REF_NAME   = config["reference_name"]    # full name e.g. "Genus_sps1"

    # ── Discover all FASTA files in the original fasta_dir ───────────────────
    _fasta_files = []
    for _ext in (".fa", ".fna", ".fasta"):
        _fasta_files.extend(_glob.glob(os.path.join(_FASTA_DIR, f"*{_ext}")))

    if not _fasta_files:
        raise ValueError(
            f"[satsuma_alignment] No FASTA files found in {_FASTA_DIR}. "
            "Check satsuma_align.fasta_dir in pipeline_paths.yaml."
        )

    _FASTA_MAP = {
        os.path.splitext(os.path.basename(f))[0]: f
        for f in _fasta_files
    }

    if _REF_NAME not in _FASTA_MAP:
        raise ValueError(
            f"[satsuma_alignment] Reference FASTA '{_REF_NAME}' not found in "
            f"{_FASTA_DIR}. Available: {list(_FASTA_MAP.keys())}"
        )

    _REF_FASTA     = _FASTA_MAP[_REF_NAME]
    _QUERY_SPECIES = sorted(k for k in _FASTA_MAP if k != _REF_NAME)

    if not _QUERY_SPECIES:
        raise ValueError(
            f"[satsuma_alignment] Only the reference FASTA was found in "
            f"{_FASTA_DIR}. At least one query species is required."
        )

    print(
        f"[satsuma_alignment] Reference    : {_REF_NAME}  ({_REF_FASTA})\n"
        f"[satsuma_alignment] Queries      : {_QUERY_SPECIES}\n"
        f"[satsuma_alignment] Threads/job  : {_THREADS}\n"
        f"[satsuma_alignment] Max parallel : {_MAX_PARALLEL}"
    )

    # ── Rule: align one query species against the reference ──────────────────
    # Output pattern uses only {query} as wildcard.
    # _REF_NAME is baked into the path string at Python parse time — that is
    # intentional and correct; Snakemake only needs to resolve {query}.
    rule run_satsuma:
        input:
            query_fa = lambda wildcards: _FASTA_MAP[wildcards.query],
            ref_fa   = _REF_FASTA
        output:
            # Single output file — the flat alignment txt the sort rule reads.
            # Named {query}.txt — just the query species name, no reference suffix.
            alignment_txt = _OUTPUT_DIR + "/{query}.txt"
        threads:
            # Snakemake reserves this many cores per job and scales down if
            # total available cores < threads * parallel_jobs.
            _THREADS
        resources:
            # satsuma_jobs is a custom token pool of size max_parallel.
            # Each job consumes 1 token — so at most max_parallel jobs
            # run simultaneously regardless of how many cores are free.
            # Pass the pool size to sybr.sh with:
            #   --resources satsuma_jobs=<max_parallel>
            # sybr.sh does this automatically from the config value.
            satsuma_jobs = 1
        params:
            tool        = _TOOL,
            ref_name    = _REF_NAME,
            out_dir     = lambda wildcards: os.path.join(
                _OUTPUT_DIR, wildcards.query
            ),
            # SATSUMA2_PATH = directory containing all Satsuma binaries.
            # Derived automatically from tool_path (same dir as the binary).
            # Override with satsuma_align.satsuma2_path in pipeline_paths.yaml.
            satsuma2_path = config["satsuma_align"].get(
                "satsuma2_path",
                os.path.dirname(config["satsuma_align"]["tool_path"])
            )
        log:
            _OUTPUT_DIR + "/logs/{query}.log"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.out_dir}
            mkdir -p $(dirname {log})

            echo "=== SatsumaSynteny2 alignment ===" | tee {log}
            echo "Query   : {input.query_fa}"        | tee -a {log}
            echo "Target  : {input.ref_fa}"          | tee -a {log}
            echo "Outdir  : {params.out_dir}"        | tee -a {log}
            echo "Threads : {threads}"               | tee -a {log}
            echo "Tool    : {params.tool}"           | tee -a {log}

            if [ ! -f "{params.tool}" ]; then
                echo "ERROR: SatsumaSynteny2 not found at {params.tool}" >&2
                echo "Check satsuma_align.tool_path in pipeline_paths.yaml"  >&2
                exit 1
            fi

            # SatsumaSynteny2 requires SATSUMA2_PATH to point to its binary dir.
            # Derived from tool_path; override with satsuma_align.satsuma2_path.
            export SATSUMA2_PATH="{params.satsuma2_path}"
            echo "SATSUMA2_PATH: $SATSUMA2_PATH" | tee -a {log}

            SATSUMA_OUT="{params.out_dir}/satsuma_summary.chained.out"

            # If the chained output already exists (e.g. from a previous interrupted
            # run), skip re-alignment and go straight to the copy step.
            if [ -f "$SATSUMA_OUT" ]; then
                echo "Resuming: satsuma_summary.chained.out already exists, skipping alignment" | tee -a {log}
            else
                "{params.tool}" \
                    -q "{input.query_fa}" \
                    -t "{input.ref_fa}"   \
                    -o "{params.out_dir}" \
                    -threads {threads} \
                    2>&1 | tee -a {log}
            fi

            if [ ! -f "$SATSUMA_OUT" ]; then
                echo "ERROR: satsuma_summary.chained.out not found in {params.out_dir}" >&2
                echo "Files present in output dir:" >&2
                ls "{params.out_dir}/" >&2
                tail -20 {log} >&2
                exit 1
            fi

            cp "$SATSUMA_OUT" "{output.alignment_txt}"

            LINE_COUNT=$(wc -l < "{output.alignment_txt}")
            echo "Alignment complete — $LINE_COUNT lines" | tee -a {log}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m Satsuma: \\e[36;1m{wildcards.query}\\e[0m vs {params.ref_name}" >&2
            """


    # ── Aggregate: wait for all species then write a marker ──────────────────
    rule satsuma_alignment_all:
        input:
            expand(
                _OUTPUT_DIR + "/{query}.txt",
                query=_QUERY_SPECIES
            )
        output:
            marker = touch(_OUTPUT_DIR + "/.satsuma_alignment_done")
        shell:
            """
            echo "All Satsuma alignments complete."
            ls -lh {input}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m all Satsuma alignments" >&2
            """
