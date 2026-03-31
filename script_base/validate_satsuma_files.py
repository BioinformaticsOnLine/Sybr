#!/usr/bin/env python3
"""
validate_satsuma_files.py — Input data verifier for the SYBR pipeline.

Runs two layers of checks:
  1. DEPENDENCY CHECK  — ensures upstream module results exist before a
                         downstream-only module is attempted.
  2. FORMAT CHECKS     — file existence, column counts, value ranges, etc.
                         (only for stages that are enabled in the config)

Dependency graph
────────────────
  synteny_processing  ──►  eba_analysis  ──►  enrichment_analysis
  chainNet_generation ──►  Ancestor_seq_recunstruction

Exit codes:
    0  — all checks passed
    1  — one or more checks failed (dependency OR format)
"""

import argparse
import os
import sys
import yaml
import glob


# ─────────────────────────────────────────────────────────────────────────────
#  ANSI colours
# ─────────────────────────────────────────────────────────────────────────────
GREEN  = "\033[0;32m"
RED    = "\033[0;31m"
YELLOW = "\033[0;33m"
CYAN   = "\033[0;36m"
BOLD   = "\033[1m"
NC     = "\033[0m"

def ok(msg):    print(f"  {GREEN}✓{NC} {msg}")
def fail(msg):  print(f"  {RED}✗{NC} {msg}", file=sys.stderr)
def warn(msg):  print(f"  {YELLOW}⚠{NC} {msg}")
def info(msg):  print(f"  {CYAN}→{NC} {msg}")


# ─────────────────────────────────────────────────────────────────────────────
#  Config helpers
# ─────────────────────────────────────────────────────────────────────────────
def load_yaml(path):
    try:
        with open(path) as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        print(f"{RED}[ERROR]{NC} Cannot read {path}: {e}", file=sys.stderr)
        sys.exit(1)


def deep_merge(base, override):
    result = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = {**result[k], **v}
        else:
            result[k] = v
    return result


def rebase_paths(cfg, in_base, out_base):
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


def build_config(merged_config_path=None, paths_file=None, user_file=None):
    if merged_config_path:
        cfg = load_yaml(merged_config_path)
    else:
        base = load_yaml(paths_file)
        user = load_yaml(user_file)
        cfg  = deep_merge(base, user)

    in_base  = cfg.get("base_input_dir",  "inputs")
    out_base = cfg.get("base_output_dir", "outputs")
    if in_base != "inputs" or out_base != "outputs":
        cfg = rebase_paths(cfg, in_base, out_base)
    return cfg


# ─────────────────────────────────────────────────────────────────────────────
#  Generic file / directory helpers
# ─────────────────────────────────────────────────────────────────────────────
def check_file(path, label, min_lines=1):
    if not os.path.isfile(path):
        fail(f"{label}: not found  →  {path}")
        return False
    size_kb = os.path.getsize(path) / 1024
    lines   = sum(1 for _ in open(path, errors="replace"))
    if lines < min_lines:
        fail(f"{label}: too few lines ({lines} < {min_lines})  →  {path}")
        return False
    ok(f"{label}  ({lines} lines, {size_kb:.1f} KB)")
    return True


def check_dir(path, label, min_files=0, ext_glob="*"):
    if not os.path.isdir(path):
        fail(f"{label}: directory not found  →  {path}")
        return False
    files = glob.glob(os.path.join(path, ext_glob))
    if len(files) < min_files:
        fail(f"{label}: only {len(files)} file(s) found (need ≥ {min_files})  →  {path}")
        return False
    ok(f"{label}  ({len(files)} file(s))")
    return True


def check_tsv_columns(path, label, expected_cols, sep="\t"):
    errs = []
    with open(path, errors="replace") as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(sep)
            if len(parts) != expected_cols:
                errs.append(f"line {i}: expected {expected_cols} cols, got {len(parts)}")
            if len(errs) >= 5:
                errs.append("… (more errors omitted)")
                break
    if errs:
        fail(f"{label}: column-count errors in {path}")
        for e in errs:
            print(f"    {e}", file=sys.stderr)
        return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
#  LAYER 1 — Dependency / prerequisite checks
# ─────────────────────────────────────────────────────────────────────────────

def _synteny_outputs_exist(cfg):
    """True if synteny_processing produced its final marker outputs."""
    synteny_results = cfg.get("synteny_results", "")
    out_final       = cfg.get("out_final", "pre-EBA1")
    if not synteny_results:
        return False
    done_files = glob.glob(
        os.path.join(synteny_results, "synteny_out", "*", "*_out", "synteny_assign_done")
    )
    pre_eba_marker = os.path.join(synteny_results, out_final, ".done_marker")
    return bool(done_files) and os.path.isfile(pre_eba_marker)


def _eba_outputs_exist(cfg):
    """True if eba_analysis produced its key output files."""
    eba_fmt = cfg.get("eba_format", {})
    pre_EBA = eba_fmt.get("pre_EBA_dir", "")
    mshsbs  = eba_fmt.get("mshsbs_dir", "")
    result  = os.path.join(pre_EBA, "EBA_WORKING", "EBA_OUT", "Merge", "Result_Merge.final") if pre_EBA else ""
    msfile  = os.path.join(mshsbs, "msHSBs.txt") if mshsbs else ""
    return (os.path.isfile(result) if result else False) and \
           (os.path.isfile(msfile) if msfile else False)


def _chainnet_outputs_exist(cfg):
    """True if chainNet_generation produced its sentinel done-files."""
    out_dir = cfg.get("chainNet", {}).get("output_dir", "")
    if not out_dir:
        return False
    sentinels = [
        os.path.join(out_dir, "clean_axt.done"),
        os.path.join(out_dir, "netSyntenic.done"),
    ]
    return all(os.path.isfile(s) for s in sentinels)


def check_dependencies(cfg, run_stages):
    """
    Enforce the two dependency chains:
      synteny_processing → eba_analysis → enrichment_analysis
      chainNet_generation → Ancestor_seq_recunstruction

    Returns the number of broken-dependency errors found.
    """
    errors = 0

    synteny_on = run_stages.get("synteny_processing",         False)
    eba_on     = run_stages.get("eba_analysis",               False)
    enrich_on  = run_stages.get("enrichment_analysis",        False)
    chain_on   = run_stages.get("chainNet_generation",         False)
    deschr_on  = run_stages.get("Ancestor_seq_recunstruction", False)

    # ── Chain A ──────────────────────────────────────────────────────────────
    # eba_analysis needs synteny_processing
    if eba_on and not synteny_on:
        if _synteny_outputs_exist(cfg):
            ok("synteny_processing outputs found from a previous run  "
               "→  eba_analysis can proceed")
        else:
            fail(
                "eba_analysis: true  but  synteny_processing: false\n"
                f"  {RED}No synteny_processing outputs found in the output directory.{NC}\n"
                f"  {YELLOW}Fix:{NC} set  synteny_processing: true  in run_sybr_config.yaml,\n"
                f"       OR point base_output_dir to a directory that already contains\n"
                f"       the synteny_processing results from a previous run."
            )
            errors += 1

    # enrichment_analysis needs eba_analysis
    if enrich_on and not eba_on:
        if _eba_outputs_exist(cfg):
            ok("eba_analysis outputs found from a previous run  "
               "→  enrichment_analysis can proceed")
        else:
            fail(
                "enrichment_analysis: true  but  eba_analysis: false\n"
                f"  {RED}No eba_analysis outputs found in the output directory.{NC}\n"
                f"  {YELLOW}Fix:{NC} set  eba_analysis: true  in run_sybr_config.yaml,\n"
                f"       OR point base_output_dir to a directory that already contains\n"
                f"       the eba_analysis results from a previous run."
            )
            errors += 1

    # enrichment also needs synteny (indirectly); surface this if eba is missing too
    if enrich_on and not eba_on and not synteny_on:
        if not _synteny_outputs_exist(cfg):
            info(
                "Note: enrichment_analysis also requires synteny_processing outputs "
                "(indirectly via eba_analysis)."
            )

    # ── Chain B ──────────────────────────────────────────────────────────────
    # Ancestor_seq_recunstruction needs chainNet_generation
    if deschr_on and not chain_on:
        if _chainnet_outputs_exist(cfg):
            ok("chainNet_generation outputs found from a previous run  "
               "→  Ancestor_seq_reconstruction can proceed")
        else:
            fail(
                "Ancestor_seq_recunstruction: true  but  chainNet_generation: false\n"
                f"  {RED}No chainNet_generation outputs found in the output directory.{NC}\n"
                f"  {YELLOW}Fix:{NC} set  chainNet_generation: true  in run_sybr_config.yaml,\n"
                f"       OR point base_output_dir to a directory that already contains\n"
                f"       the chainNet_generation results from a previous run."
            )
            errors += 1

    if errors == 0 and not (
        (eba_on and not synteny_on) or
        (enrich_on and not eba_on)  or
        (deschr_on and not chain_on)
    ):
        ok("All module dependencies satisfied")

    return errors


# ─────────────────────────────────────────────────────────────────────────────
#  LAYER 2 — Module-specific format validators
# ─────────────────────────────────────────────────────────────────────────────

def validate_satsuma(cfg, verbose):
    errors      = 0
    satsuma_dir = cfg.get("satsuma_alignments", "")
    if not os.path.isdir(satsuma_dir):
        fail(f"Satsuma dir not found: {satsuma_dir}")
        return 1

    txt_files = sorted(glob.glob(os.path.join(satsuma_dir, "*.txt")))
    if not txt_files:
        fail(f"No .txt files found in {satsuma_dir}")
        return 1

    print(f"\n  {'='*58}")
    print(f"  {'SATSUMA ALIGNMENT FILES':^58}")
    print(f"  {'='*58}")

    total_lines, total_kb = 0, 0.0
    for fpath in txt_files:
        fname   = os.path.basename(fpath)
        lines   = sum(1 for _ in open(fpath, errors="replace"))
        size_kb = os.path.getsize(fpath) / 1024
        total_lines += lines
        total_kb    += size_kb

        sample_errs = []
        with open(fpath, errors="replace") as f:
            for i, line in enumerate(f, 1):
                if i > 200:
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 8:
                    sample_errs.append(f"line {i}: {len(parts)} cols (expected 8)")
                    continue
                try:
                    float(parts[6])
                except ValueError:
                    sample_errs.append(f"line {i}: col 7 not numeric ('{parts[6]}')")
                if parts[7] not in ("+", "-"):
                    sample_errs.append(f"line {i}: col 8 not +/- ('{parts[7]}')")
                if len(sample_errs) >= 3:
                    break

        sym = f"{GREEN}✓{NC}" if not sample_errs else f"{RED}✗{NC}"
        print(f"  {sym} {fname:<35} {lines:>8} lines {size_kb:>10.1f} KB")
        if sample_errs and verbose:
            for e in sample_errs[:3]:
                print(f"       {YELLOW}→ {e}{NC}")
        if sample_errs:
            errors += 1

    print(f"  {'-'*58}")
    print(f"  Total: {len(txt_files):>3} files  {total_lines:>10} lines  {total_kb:>10.1f} KB")
    print(f"  {'='*58}")
    return errors


def validate_sequence_lengths(cfg):
    path = cfg.get("sequence_lengths_file", "")
    if not path:
        warn("sequence_lengths_file not set — skipping")
        return 0
    if not check_file(path, "all_sequence_lengths.txt", min_lines=1):
        return 1
    if not check_tsv_columns(path, "all_sequence_lengths.txt", expected_cols=3):
        return 1
    errs = []
    with open(path, errors="replace") as f:
        for i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                int(parts[1])
            except ValueError:
                errs.append(f"line {i}: col 2 not integer ('{parts[1]}')")
            if len(errs) >= 3:
                break
    if errs:
        for e in errs:
            fail(f"all_sequence_lengths.txt: {e}")
        return 1
    return 0


def validate_eba(cfg):
    errors  = 0
    eba_fmt = cfg.get("eba_format", {})

    eba_input_dir = cfg.get("eba", {}).get("d", eba_fmt.get("eba_input_dir", ""))
    if eba_input_dir:
        if os.path.isdir(eba_input_dir):
            check_dir(eba_input_dir, "EBA input dir")
        # else: silently skip — EBA-input is generated by the pipeline
    else:
        warn("EBA input dir not configured")

    chr_size = eba_fmt.get("chr_size_file", "")
    if chr_size:
        if check_file(chr_size, "EBA/chr_size.txt", min_lines=1):
            if not check_tsv_columns(chr_size, "EBA/chr_size.txt", expected_cols=2):
                errors += 1
            else:
                bad = []
                with open(chr_size, errors="replace") as f:
                    for i, line in enumerate(f, 1):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        try:
                            int(parts[1])
                        except ValueError:
                            bad.append(f"line {i}: '{parts[1]}' not an integer")
                if bad:
                    for b in bad:
                        fail(f"EBA/chr_size.txt: {b}")
                    errors += 1
        else:
            errors += 1
    else:
        warn("EBA chr_size_file not configured")

    cls_file = cfg.get("eba", {}).get("c", "")
    if cls_file:
        if check_file(cls_file, "EBA/classification.eba", min_lines=1):
            has_lineage = False
            fmt_errs    = []
            with open(cls_file, errors="replace") as f:
                for i, line in enumerate(f, 1):
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if "=" not in line:
                        fmt_errs.append(f"line {i}: no '=' separator")
                        continue
                    key, _, _ = line.partition("=")
                    if key == "lineage":
                        has_lineage = True
            if not has_lineage:
                fail("EBA/classification.eba: missing 'lineage=' line")
                errors += 1
            for e in fmt_errs[:3]:
                fail(f"EBA/classification.eba: {e}")
            if fmt_errs:
                errors += 1
        else:
            errors += 1
    else:
        warn("EBA classification file not configured")

    scaffolds = eba_fmt.get("scaffolds_file", "ALL_CHROMOSOMES")
    if scaffolds and scaffolds not in ("ALL_CHROMOSOMES", "Chromosomes"):
        if not check_file(scaffolds, "Scaffolds.txt", min_lines=1):
            errors += 1

    return errors


def validate_enrichment(cfg):
    errors  = 0
    enr_cfg = cfg.get("enrichment", {})

    ann_file = enr_cfg.get("annotation_file", "")
    if ann_file:
        if check_file(ann_file, "protein_annotation.tsv", min_lines=1):
            if not check_tsv_columns(ann_file, "protein_annotation.tsv", expected_cols=5):
                errors += 1
    else:
        warn("enrichment annotation_file not configured")

    kegg_file = enr_cfg.get("kegg_file", "")
    if kegg_file:
        if check_file(kegg_file, "3kegg_annotationTOgenes.txt", min_lines=2):
            with open(kegg_file, errors="replace") as f:
                header = f.readline().rstrip("\n").split("\t")
            if header != ["term", "gene"]:
                fail(f"3kegg_annotationTOgenes.txt: unexpected header {header} (want ['term','gene'])")
                errors += 1
            elif not check_tsv_columns(kegg_file, "3kegg_annotationTOgenes.txt", expected_cols=2):
                errors += 1
    else:
        warn("enrichment kegg_file not configured")

    return errors


def validate_chainnet(cfg):
    errors    = 0
    chain_cfg = cfg.get("chainNet", {})

    lastz_dir = chain_cfg.get("lastZ_alignments", "")
    if lastz_dir:
        axt_files = glob.glob(os.path.join(lastz_dir, "*.axt"))
        if not axt_files:
            fail(f"LastZ_alignments: no .axt files in {lastz_dir}")
            errors += 1
        else:
            ok(f"LastZ_alignments  ({len(axt_files)} .axt files)")
            for axt in axt_files[:3]:
                with open(axt, errors="replace") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        if not line.split()[0].isdigit():
                            fail(f"{os.path.basename(axt)}: first data line missing block index")
                            errors += 1
                        break
    else:
        warn("chainNet lastZ_alignments not configured")

    seq_dir = chain_cfg.get("seq_dir", "")
    if seq_dir:
        fa_files = (glob.glob(os.path.join(seq_dir, "*.fa"))  +
                    glob.glob(os.path.join(seq_dir, "*.fasta")) +
                    glob.glob(os.path.join(seq_dir, "*.fna")))
        if not fa_files:
            fail(f"seq/: no FASTA files in {seq_dir}")
            errors += 1
        else:
            ok(f"seq/  ({len(fa_files)} FASTA files)")
            for fa in fa_files:
                with open(fa, errors="replace") as f:
                    first = f.readline().strip()
                if not first.startswith(">"):
                    fail(f"{os.path.basename(fa)}: does not start with '>' (invalid FASTA)")
                    errors += 1
    else:
        warn("chainNet seq_dir not configured")

    return errors


def validate_deschrambler(cfg):
    errors = 0
    dsc    = cfg.get("deschrambler", {})

    tree_file = dsc.get("tree_file", "")
    if tree_file:
        if check_file(tree_file, "tree.txt", min_lines=1):
            with open(tree_file, errors="replace") as f:
                content = f.read().strip()
            if not (content.startswith("(") and content.endswith(";")):
                fail("tree.txt: not a valid Newick tree (must start with '(' and end with ';')")
                errors += 1
    else:
        warn("deschrambler tree_file not configured")

    sp_file = dsc.get("species_file", "")
    if sp_file:
        if check_file(sp_file, "species_info.txt", min_lines=1):
            bad = []
            with open(sp_file, errors="replace") as f:
                for i, line in enumerate(f, 1):
                    if len(line.rstrip("\n").split()) < 2:
                        bad.append(f"line {i}: fewer than 2 fields")
            for b in bad[:3]:
                fail(f"species_info.txt: {b}")
            if bad:
                errors += 1
    else:
        warn("deschrambler species_file not configured")

    tool_dir = dsc.get("pipeline_tool_dir", "")
    if tool_dir:
        if not os.path.isdir(tool_dir):
            fail(f"DESCHRAMBLER tool dir not found: {tool_dir}")
            errors += 1
        else:
            ok("DESCHRAMBLER tool dir found")
    else:
        warn("deschrambler pipeline_tool_dir not configured")

    return errors


# ─────────────────────────────────────────────────────────────────────────────
#  Section header helpers
# ─────────────────────────────────────────────────────────────────────────────
def section(title):
    bar = "─" * max(0, 55 - len(title))
    print(f"\n{BOLD}{CYAN}┌─ {title} {bar}┐{NC}")

def section_result(errors):
    label = f"{GREEN}PASSED{NC}" if errors == 0 else f"{RED}FAILED ({errors} error(s)){NC}"
    print(f"{BOLD}{CYAN}└─ Result: {label}{NC}")


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="SYBR pipeline input validator")
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--config",      help="Pre-merged config YAML (from sybr.sh)")
    grp.add_argument("--paths",       help="pipeline_paths.yaml  (use with --user-config)")
    parser.add_argument("--user-config", help="run_sybr_config.yaml")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show extra detail on failures")
    args = parser.parse_args()

    if args.paths and not args.user_config:
        parser.error("--paths requires --user-config")

    cfg = (build_config(merged_config_path=args.config)
           if args.config
           else build_config(paths_file=args.paths, user_file=args.user_config))

    run_stages   = cfg.get("run_stages", {})
    total_errors = 0

    print(f"\n{BOLD}╔══════════════════════════════════════════════════════════╗")
    print(  f"║           SYBR Pipeline — Input Verification             ║")
    print(  f"╚══════════════════════════════════════════════════════════╝{NC}")

    # ── LAYER 1: dependency checks (always runs, blocks on failure) ───────────
    section("Checking module dependencies")
    dep_errors = check_dependencies(cfg, run_stages)
    section_result(dep_errors)
    total_errors += dep_errors

    if dep_errors:
        print(f"\n{RED}{BOLD}Pipeline aborted: resolve the dependency issues above "
              f"before re-running.{NC}")
        print(f"\n{'═'*60}")
        print(f"{RED}✗ Dependency check failed. Pipeline cannot run as configured.{NC}")
        print(f"{'═'*60}\n")
        sys.exit(1)

    # ── LAYER 2: per-module input format checks ───────────────────────────────
    if run_stages.get("synteny_processing", False):
        section("synteny_processing — input format checks")
        errs  = validate_satsuma(cfg, args.verbose)
        errs += validate_sequence_lengths(cfg)
        total_errors += errs
        section_result(errs)

    if run_stages.get("eba_analysis", False):
        section("eba_analysis — input format checks")
        errs = validate_eba(cfg)
        total_errors += errs
        section_result(errs)

    if run_stages.get("enrichment_analysis", False):
        section("enrichment_analysis — input format checks")
        errs = validate_enrichment(cfg)
        total_errors += errs
        section_result(errs)

    if run_stages.get("chainNet_generation", False):
        section("chainNet_generation — input format checks")
        errs = validate_chainnet(cfg)
        total_errors += errs
        section_result(errs)

    if run_stages.get("Ancestor_seq_recunstruction", False):
        section("Ancestor_seq_reconstruction — input format checks")
        errs = validate_deschrambler(cfg)
        total_errors += errs
        section_result(errs)

    # ── Grand summary ─────────────────────────────────────────────────────────
    print(f"\n{'═'*60}")
    if total_errors == 0:
        print(f"{GREEN}✓ All checks passed — pipeline is ready to run.{NC}")
    else:
        print(f"{RED}✗ {total_errors} check(s) failed. "
              f"Fix the issues above before running the pipeline.{NC}")
    print(f"{'═'*60}\n")

    sys.exit(0 if total_errors == 0 else 1)


if __name__ == "__main__":
    main()
