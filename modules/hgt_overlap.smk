# hgt_overlap.smk — HGT × msHSBs / EBRs / ancestoral_genes overlap analysis
#
# When hgt_overlap_analysis: true, overlaps hgt.txt against up to three
# reference sets depending on which upstream stages are also enabled:
#   msHSBs           — requires eba_analysis
#   EBRs             — requires eba_analysis
#   ancestoral_genes — requires Ancestor_seq_recunstruction (deschrambler block_list)
#
# Chromosome normalisation: hgt.txt uses "chr1" style; all three reference
# files use bare numbers ("1"). The "chr" prefix is stripped from hgt.txt.
# Overlap criterion: same chromosome AND coordinate ranges share ≥ 1 bp.
#
# Outputs (all in <hgt_overlap_dir>/):
#   hgt_overlapping_msHSBs.txt           / hgt_msHSBs_overlap_summary.txt
#   hgt_overlapping_EBRs.txt             / hgt_EBRs_overlap_summary.txt
#   hgt_overlapping_ancestoral_genes.txt / hgt_ancestoral_genes_overlap_summary.txt

# ── Shared helper used by every run: block ───────────────────────────────────
def _hgt_overlap_run(hgt_file, ref_file, ref_cols, out_overlapping, out_summary,
                     out_dir, ref_label):
    """
    Core overlap logic shared by all three HGT overlap rules.

    ref_cols  : number of columns to keep from each reference line (e.g. 3 for
                msHSBs/EBRs, 5 for block_list).
    ref_label : human-readable label used in the summary file.
    """
    import os
    os.makedirs(out_dir, exist_ok=True)

    # Parse reference file — strip "chr" prefix if present for consistency
    regions = []
    with open(ref_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            r_chr = p[0][3:] if p[0].lower().startswith("chr") else p[0]
            regions.append((r_chr, int(p[1]), int(p[2]), p[:ref_cols]))

    # Parse hgt.txt — strip "chr" prefix to match reference format
    overlaps      = []
    hgt_total     = 0
    hgt_overlap_n = 0

    with open(hgt_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            hgt_total += 1
            raw_chr   = parts[0]
            hgt_chr   = raw_chr[3:] if raw_chr.lower().startswith("chr") else raw_chr
            hgt_start = int(parts[1])
            hgt_end   = int(parts[2])

            for (r_chr, r_start, r_end, r_parts) in regions:
                if hgt_chr != r_chr:
                    continue
                if hgt_end < r_start or hgt_start > r_end:
                    continue
                overlap_start = max(hgt_start, r_start)
                overlap_end   = min(hgt_end,   r_end)
                overlap_len   = overlap_end - overlap_start
                overlaps.append(parts + r_parts +
                                 [str(overlap_start), str(overlap_end), str(overlap_len)])
                hgt_overlap_n += 1

    # Write overlap hits (trim/pad to actual column count)
    n_ref_cols = ref_cols          # columns from reference kept in output
    n_hgt_cols = 4                 # chr, start, end, gene (hgt.txt always has 4)
    total_cols  = n_hgt_cols + n_ref_cols + 3   # + overlap_start/end/len
    with open(out_overlapping, "w") as out:
        # Build header dynamically
        hgt_hdr = ["hgt_chr", "hgt_start", "hgt_end", "hgt_gene"]
        if ref_label == "msHSBs":
            ref_hdr = ["msHSB_chr", "msHSB_start", "msHSB_end"]
        elif ref_label == "EBRs":
            ref_hdr = ["EBR_chr", "EBR_start", "EBR_end", "EBR_species"]
        else:  # ancestoral_genes
            ref_hdr = ["block_chr", "block_start", "block_end", "block_strand", "block_id"]
        hdr = hgt_hdr + ref_hdr + ["overlap_start", "overlap_end", "overlap_length_bp"]
        out.write("\t".join(hdr) + "\n")
        for row in overlaps:
            out.write("\t".join(row[:len(hdr)]) + "\n")

    # Write summary
    with open(out_summary, "w") as out:
        out.write(f"# HGT × {ref_label} overlap summary\n")
        out.write(f"HGT file                 : {hgt_file}\n")
        out.write(f"{ref_label} file         : {ref_file}\n")
        out.write(f"Total {ref_label} regions: {len(regions)}\n")
        out.write(f"Total HGT entries        : {hgt_total}\n")
        out.write(f"HGT entries with overlap : {hgt_overlap_n}\n")
        out.write(f"Total overlap rows       : {len(overlaps)}\n")
        out.write(f"Output file              : {out_overlapping}\n")

    print(
        f"\033[32;1m✔\033[0m HGT × {ref_label} overlap complete — "
        f"{hgt_overlap_n}/{hgt_total} HGT entries overlap with {ref_label}",
        file=__import__("sys").stderr
    )


# ── Common paths (defined once, used by all rules below) ─────────────────────
if should_run_rule("hgt_overlap_analysis"):

    _HGT_OVERLAP_INPUT_DIR = config.get("hgt_overlap", {}).get(
        "hgt_dir", os.path.join(config.get("base_input_dir", "inputs"), "HGTs")
    )
    _HGT_FILE = os.path.join(_HGT_OVERLAP_INPUT_DIR, "hgt.txt")

    _HGT_OVERLAP_OUT_DIR = config.get("hgt_overlap", {}).get(
        "output_dir",
        os.path.join(config.get("base_output_dir", "outputs"), "hgt_overlap")
    )

    # ── Rule 1: HGT × msHSBs  (needs eba_analysis) ──────────────────────────
    if should_run_rule("eba_analysis"):

        _MSHSBS_DIR  = config["eba_format"]["mshsbs_dir"]
        _MSHSBS_FILE = os.path.join(_MSHSBS_DIR, "msHSBs.txt")
        _MSHSBS_DONE = os.path.join(_MSHSBS_DIR, "msHSBs_done.txt")

        rule hgt_overlap_msHSBs:
            input:
                msHSBs_done = _MSHSBS_DONE,
                hgt_file    = _HGT_FILE,
            output:
                overlapping = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_overlapping_msHSBs.txt"),
                summary     = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_msHSBs_overlap_summary.txt"),
            params:
                out_dir     = _HGT_OVERLAP_OUT_DIR,
                ref_file    = _MSHSBS_FILE,
            resources:
                mem_mb=4000, time_min=15
            run:
                import sys as _sys
                print("\n\033[36;1m━━━  HGT Overlap Analysis  ━━━\033[0m", file=_sys.stderr)
                print("\033[36m▶  Overlapping HGT coordinates with msHSBs …\033[0m", file=_sys.stderr)
                _hgt_overlap_run(
                    hgt_file      = input.hgt_file,
                    ref_file      = params.ref_file,
                    ref_cols      = 3,
                    out_overlapping = output.overlapping,
                    out_summary   = output.summary,
                    out_dir       = params.out_dir,
                    ref_label     = "msHSBs",
                )

    # ── Rule 2: HGT × EBRs  (needs eba_analysis) ────────────────────────────
    if should_run_rule("eba_analysis"):

        _EBRS_DIR  = config["eba_format"]["ebrs_dir"]
        _EBRS_FILE = os.path.join(_EBRS_DIR, "EBRs.txt")
        _EBRS_DONE = os.path.join(_EBRS_DIR, "ebrs_processed.marker")

        rule hgt_overlap_EBRs:
            input:
                ebrs_done = _EBRS_DONE,
                hgt_file  = _HGT_FILE,
            output:
                overlapping = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_overlapping_EBRs.txt"),
                summary     = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_EBRs_overlap_summary.txt"),
            params:
                out_dir  = _HGT_OVERLAP_OUT_DIR,
                ref_file = _EBRS_FILE,
            resources:
                mem_mb=4000, time_min=15
            run:
                import sys as _sys
                print("\033[36m▶  Overlapping HGT coordinates with EBRs …\033[0m", file=_sys.stderr)
                _hgt_overlap_run(
                    hgt_file      = input.hgt_file,
                    ref_file      = params.ref_file,
                    ref_cols      = 4,
                    out_overlapping = output.overlapping,
                    out_summary   = output.summary,
                    out_dir       = params.out_dir,
                    ref_label     = "EBRs",
                )

    # ── Rule 3: HGT × ancestoral_genes  (needs Ancestor_seq_recunstruction) ─────
    if should_run_rule("Ancestor_seq_recunstruction"):

        _ALIGN_OUTPUT = config["chainNet"].get("output_dir",
                            "outputs/Ancestor_seq_recunstruction")
        _BLOCK_LIST   = os.path.join(_ALIGN_OUTPUT, "APCFs.300K", "SFs", "block_list.txt")
        _DESCHRAMBLER_DONE = f"{_ALIGN_OUTPUT}/deschrambler.done"

        rule hgt_overlap_ancestoral_genes:
            input:
                deschrambler_done = _DESCHRAMBLER_DONE,
                hgt_file          = _HGT_FILE,
            output:
                overlapping = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_overlapping_ancestoral_genes.txt"),
                summary     = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_ancestoral_genes_overlap_summary.txt"),
            params:
                out_dir  = _HGT_OVERLAP_OUT_DIR,
                ref_file = _BLOCK_LIST,
            resources:
                mem_mb=4000, time_min=15
            run:
                import sys as _sys
                print("\033[36m▶  Overlapping HGT coordinates with ancestoral genes (DESCHRAMBLER block list) …\033[0m", file=_sys.stderr)
                _hgt_overlap_run(
                    hgt_file      = input.hgt_file,
                    ref_file      = params.ref_file,
                    ref_cols      = 5,
                    out_overlapping = output.overlapping,
                    out_summary   = output.summary,
                    out_dir       = params.out_dir,
                    ref_label     = "ancestoral_genes",
                )

    # ── Rule 4: Interactive HTML tree (needs both EBRs + ancestoral_genes rules) ─
    if should_run_rule("eba_analysis") and should_run_rule("Ancestor_seq_recunstruction"):

        _TREE_FILE   = config.get("deschrambler", {}).get(
            "tree_file",
            os.path.join(config.get("base_input_dir", "inputs"),
                         "Ancestor_seq_recunstruction", "tree.txt")
        )
        _SCRIPTS_DIR = config.get("scripts", "script_base")
        _TREE_SCRIPT = os.path.join(_SCRIPTS_DIR, "generate_hgt_tree.py")

        rule hgt_tree_html:
            input:
                tree           = _TREE_FILE,
                ebr_overlaps   = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_overlapping_EBRs.txt"),
                block_overlaps = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_overlapping_ancestoral_genes.txt"),
            output:
                html = os.path.join(_HGT_OVERLAP_OUT_DIR, "hgt_tree.html"),
            params:
                script = _TREE_SCRIPT,
            resources:
                mem_mb=1000, time_min=5
            shell:
                """
                printf '\033[36m▶  Generating interactive HGT phylogenetic tree HTML …\033[0m\n' >&2
                python3 {params.script} \
                    --tree           {input.tree} \
                    --ebr-overlaps   {input.ebr_overlaps} \
                    --block-overlaps {input.block_overlaps} \
                    --output         {output.html}
                printf '\033[32;1m✔\033[0m  Interactive HTML tree ready → %s\n' '{output.html}' >&2
                printf '\033[33m    Features: 5 layout modes · zoom/pan · light/dark theme · gene search\033[0m\n' >&2
                """
