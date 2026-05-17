# ─────────────────────────────────────────────────────────────────────────────
#  lastz_alignment.smk  —  Optional de-novo LastZ alignment module
#
#  Mirrors the structure of satsuma_alignment.smk exactly.
#
#  Activated by setting  run_stages.run_lastz_alignment: true  in
#  run_sybr_config.yaml.
#
#  When active:
#    - FASTA files are read from  satsuma_align.fasta_dir  (same source as
#      satsuma_alignment.smk — one FASTA per species, full names as stems)
#    - The reference is identified by extracting the species part of
#      reference_name  (second word, lowercased):
#        e.g.  "Adineta_vaga"  →  "vaga"
#      This matches the short-name convention used for .axt file stems
#        e.g.  habrotrocha.axt  ricciae.axt  roseola.axt
#    - LastZ is run for every non-reference species, producing {species}.axt
#      directly in  chainNet.lastZ_alignments
#    - alignment_processing.smk and deschrambler.smk then pick up those
#      .axt files as normal
#
#  When inactive:
#    - This module defines no rules and does nothing
#    - alignment_processing.smk reads pre-computed .axt files as before
#
#  LastZ command:
#    lastz reference.fa[multiple] query.fa \
#        C=0 E=30 H=2000 K=2200 L=2200 O=400 Y=3400 \
#        --format=axt --output={species}.axt
# ─────────────────────────────────────────────────────────────────────────────

import os
import glob as _glob

if should_run_rule("run_lastz_alignment"):

    # ── Resolve paths and settings from config ───────────────────────────────
    # chr-FASTA files use short names (sps1.fa, not Genus_sps1.fa) and
    # >chr* headers, which LastZ requires.
    _FASTA_DIR    = CHR_FASTA_DIR   # converted copies live here
    # config_handling.smk has already overridden chainNet.lastZ_alignments to
    # point at lastz_align.output_dir (an outputs/ path) when
    # run_lastz_alignment is true, so this always resolves to the correct dir.
    _AXT_DIR      = config["chainNet"]["lastZ_alignments"]
    _LASTZ_BIN    = config["lastz_align"]["tool_path"]
    _LASTZ_PARAMS = config["lastz_align"].get(
        "params", "C=0 E=30 H=2000 K=2200 L=2200 O=400 Y=3400"
    )
    _MAX_PARALLEL = config["lastz_align"].get("max_parallel", 2)
    _THREADS      = config["lastz_align"].get("threads", 8)

    # ── Derive the short reference name ─────────────────────────────────────
    # In chr_fasta_dir filenames are already short: sps1.fa (no Genus_ prefix).
    # reference_species gives the short name directly (e.g. "sps1").
    _REF_FULL  = config["reference_name"]           # e.g. "Genus_sps1"
    _REF_SHORT = config["reference_species"]        # e.g. "sps1" — matches chr_fasta_dir filename

    # ── Discover species from the ORIGINAL fasta_dir (always present) ───────
    # CHR_FASTA_DIR doesn't exist yet at parse time — prepare_chr_fasta creates
    # it at runtime.  We therefore derive species names and expected chr-FASTA
    # paths from the original fasta_dir (which always exists), exactly as
    # alignment_processing.smk does.
    _ORIG_FASTA_DIR_LZ = config["satsuma_align"]["fasta_dir"]
    _orig_fasta_files = []
    for _ext in (".fa", ".fna", ".fasta"):
        _orig_fasta_files.extend(_glob.glob(os.path.join(_ORIG_FASTA_DIR_LZ, f"*{_ext}")))

    if not _orig_fasta_files:
        raise ValueError(
            f"[lastz_alignment] No FASTA files found in original fasta_dir: "
            f"{_ORIG_FASTA_DIR_LZ}. Check satsuma_align.fasta_dir in pipeline_paths.yaml."
        )

    # Build short_stem → expected chr-FASTA path map.
    # Convention: short stem = everything after first underscore (Genus_sps1 → sps1).
    # If no underscore, use the full stem as-is.
    _FASTA_MAP_FULL = {}
    for _f in _orig_fasta_files:
        _full_stem = os.path.splitext(os.path.basename(_f))[0]   # e.g. "Genus_sps1"
        _short     = _full_stem.split("_", 1)[-1]                # e.g. "sps1"
        _chr_path  = os.path.join(CHR_FASTA_DIR, _short + ".fa") # expected chr-FASTA
        _FASTA_MAP_FULL[_short] = _chr_path

    # Verify the reference short name maps to a known species
    if _REF_SHORT not in _FASTA_MAP_FULL:
        raise ValueError(
            f"[lastz_alignment] Reference species '{_REF_SHORT}' not found among "
            f"original FASTA stems in {_ORIG_FASTA_DIR_LZ}. "
            f"Available short names: {list(_FASTA_MAP_FULL.keys())}"
        )

    _REF_FASTA = _FASTA_MAP_FULL[_REF_SHORT]   # expected path in CHR_FASTA_DIR

    # ── Build query map ──────────────────────────────────────────────────────
    _QUERY_MAP = {
        stem: path
        for stem, path in _FASTA_MAP_FULL.items()
        if stem != _REF_SHORT
    }

    if not _QUERY_MAP:
        raise ValueError(
            f"[lastz_alignment] No query species found in {_ORIG_FASTA_DIR_LZ} "
            f"(excluding reference '{_REF_FULL}')."
        )

    _QUERY_SPECIES = sorted(_QUERY_MAP.keys())

    print(
        f"[lastz_alignment] Reference      : {_REF_FULL} → short name '{_REF_SHORT}'\n"
        f"[lastz_alignment] Reference FASTA: {_REF_FASTA}\n"
        f"[lastz_alignment] Query species  : {_QUERY_SPECIES}\n"
        f"[lastz_alignment] AXT output dir : {_AXT_DIR}\n"
        f"[lastz_alignment] Threads/job    : {_THREADS}\n"
        f"[lastz_alignment] Max parallel   : {_MAX_PARALLEL}\n"
        f"[lastz_alignment] Params         : {_LASTZ_PARAMS}"
    )

    # Ensure AXT output directory exists
    os.makedirs(_AXT_DIR, exist_ok=True)

    # ── Rule: run LastZ for one query species ────────────────────────────────
    rule run_lastz:
        input:
            query_fa  = lambda wildcards: _QUERY_MAP[wildcards.species],
            ref_fa    = _REF_FASTA,
            chr_ready = CHR_FASTA_DIR + "/.chr_fasta_done"
        output:
            axt = _AXT_DIR + "/{species}.axt"
        threads:
            _THREADS
        resources:
            lastz_jobs = 1
        params:
            tool      = _LASTZ_BIN,
            lz_params = _LASTZ_PARAMS,
            axt_dir   = _AXT_DIR,
            ref_name  = _REF_FULL
        log:
            _AXT_DIR + "/logs/{species}_lastz.log"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.axt_dir}
            mkdir -p $(dirname {log})

            echo "=== LastZ alignment ===" | tee {log}
            echo "Query     : {input.query_fa}"   | tee -a {log}
            echo "Reference : {input.ref_fa}"     | tee -a {log}
            echo "Output    : {output.axt}"       | tee -a {log}
            echo "Threads   : {threads}"          | tee -a {log}
            echo "Tool      : {params.tool}"      | tee -a {log}
            echo "Params    : {params.lz_params}" | tee -a {log}

            # Resolve the lastz binary — use configured path if it is a file,
            # otherwise fall back to lastz on PATH (e.g. conda install)
            LASTZ_CMD="{params.tool}"
            if [ ! -f "$LASTZ_CMD" ]; then
                if command -v lastz &>/dev/null; then
                    LASTZ_CMD="lastz"
                else
                    echo "ERROR: lastz not found at '{params.tool}' and not on PATH" >&2
                    echo "Check lastz_align.tool_path in pipeline_paths.yaml" >&2
                    exit 1
                fi
            fi

            "$LASTZ_CMD" \
                "{input.ref_fa}[multiple]" \
                "{input.query_fa}" \
                {params.lz_params} \
                --format=axt \
                --output="{output.axt}" \
                2>&1 | tee -a {log}

            if [ ! -s "{output.axt}" ]; then
                echo "ERROR: LastZ produced an empty or missing .axt file" >&2
                exit 1
            fi

            LINE_COUNT=$(wc -l < "{output.axt}")
            echo "Alignment complete — $LINE_COUNT lines in .axt" | tee -a {log}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m LastZ: \\e[36;1m{wildcards.species}\\e[0m vs {params.ref_name}" >&2
            """


    # ── Aggregate: wait for all .axt files then write a marker ───────────────
    rule lastz_alignment_all:
        input:
            expand(_AXT_DIR + "/{species}.axt", species=_QUERY_SPECIES)
        output:
            marker = touch(_AXT_DIR + "/.lastz_alignment_done")
        shell:
            """
            echo "All LastZ alignments complete."
            ls -lh {input}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m all LastZ alignments" >&2
            """
