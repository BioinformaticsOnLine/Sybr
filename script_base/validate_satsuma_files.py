#!/usr/bin/env python3
"""
validate_satsuma_files.py
=========================
Module-aware input validator for the SYBR pipeline.

Checks both existence AND data format of required input files.

Required inputs per module combination
───────────────────────────────────────
fasta/ dir:
  any of: run_satsuma_alignment, run_lastz_alignment,
          chainNet_generation, Ancestor_seq_recunstruction  ON

synteny_processing/Satsuma_alignments/ (.txt):
  run_satsuma_alignment: false  AND any of synteny_processing,
  eba_analysis, enrichment_analysis  ON

synteny_processing/all_sequence_lengths.txt:
  synteny_processing: true

eba_analysis/classification.eba:
  eba_analysis: true

enrichment_analysis/protein_annotation.tsv:
  enrichment_analysis: true

Ancestor_seq_recunstruction/LastZ_alignments/ (.axt):
  run_lastz_alignment: false  AND any of chainNet_generation,
  Ancestor_seq_recunstruction  ON

Ancestor_seq_recunstruction/tree.txt:
  Ancestor_seq_recunstruction: true

Ancestor_seq_recunstruction/species_info.txt:
  Ancestor_seq_recunstruction: true  OR  chainNet_generation: true
"""

import argparse, glob, os, re, sys

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

GREEN  = "\033[0;32m"; RED    = "\033[0;31m"; CYAN   = "\033[0;36m"
YELLOW = "\033[0;33m"; BOLD   = "\033[1m";    NC     = "\033[0m"
CHECK  = f"{GREEN}✔{NC}"; CROSS = f"{RED}✗{NC}"; WARN = f"{YELLOW}⚠{NC}"

# Maximum lines sampled for format checking (avoids reading huge files fully)
FORMAT_SAMPLE_LINES = 50

# Valid nucleotide characters in FASTA sequences (IUPAC + gaps)
_VALID_NT = set("ACGTacgtNnRrYySsWwKkMmBbDdHhVv-\n\r")


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def has_ext(directory, *exts):
    for ext in exts:
        if glob.glob(os.path.join(directory, f"*{ext}")):
            return True
    return False


# ── Format validators ────────────────────────────────────────────────────────

def _check_satsuma_format(path):
    """
    Satsuma alignment .txt — 8 tab-separated columns:
      col0  ref_chr      integer
      col1  ref_start    integer >= 0
      col2  ref_end      integer > ref_start
      col3  query_chr    integer
      col4  query_start  integer >= 0
      col5  query_end    integer > query_start
      col6  identity     float in [0.0, 1.0]
      col7  strand       '+' or '-'
    """
    errors = []
    n = 0
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                if lineno > FORMAT_SAMPLE_LINES:
                    break
                line = raw.rstrip("\n\r")
                if not line:
                    continue
                n += 1
                cols = line.split("\t")
                if len(cols) != 8:
                    errors.append(f"  line {lineno}: expected 8 tab-separated columns, got {len(cols)}")
                    if len(errors) >= 5:
                        errors.append("  ... (further column errors suppressed)")
                        break
                    continue
                # col0 ref_chr
                try:
                    int(cols[0])
                except ValueError:
                    errors.append(f"  line {lineno}: col1 (ref_chr) must be integer, got '{cols[0]}'")
                # col1 ref_start >= 0
                try:
                    rs = int(cols[1])
                    if rs < 0:
                        errors.append(f"  line {lineno}: col2 (ref_start) must be >= 0, got {rs}")
                except ValueError:
                    errors.append(f"  line {lineno}: col2 (ref_start) must be integer, got '{cols[1]}'")
                # col2 ref_end > ref_start
                try:
                    re_ = int(cols[2]); rs_ = int(cols[1])
                    if re_ <= rs_:
                        errors.append(f"  line {lineno}: col3 (ref_end={re_}) must be > col2 (ref_start={rs_})")
                except ValueError:
                    errors.append(f"  line {lineno}: col3 (ref_end) must be integer, got '{cols[2]}'")
                # col3 query_chr
                try:
                    int(cols[3])
                except ValueError:
                    errors.append(f"  line {lineno}: col4 (query_chr) must be integer, got '{cols[3]}'")
                # col4 query_start >= 0
                try:
                    qs = int(cols[4])
                    if qs < 0:
                        errors.append(f"  line {lineno}: col5 (query_start) must be >= 0, got {qs}")
                except ValueError:
                    errors.append(f"  line {lineno}: col5 (query_start) must be integer, got '{cols[4]}'")
                # col5 query_end > query_start
                try:
                    qe = int(cols[5]); qs_ = int(cols[4])
                    if qe <= qs_:
                        errors.append(f"  line {lineno}: col6 (query_end={qe}) must be > col5 (query_start={qs_})")
                except ValueError:
                    errors.append(f"  line {lineno}: col6 (query_end) must be integer, got '{cols[5]}'")
                # col6 identity in [0,1]
                try:
                    ident = float(cols[6])
                    if not (0.0 <= ident <= 1.0):
                        errors.append(f"  line {lineno}: col7 (identity={ident}) must be in [0.0, 1.0]")
                except ValueError:
                    errors.append(f"  line {lineno}: col7 (identity) must be a float, got '{cols[6]}'")
                # col7 strand
                if cols[7] not in ("+", "-"):
                    errors.append(f"  line {lineno}: col8 (strand) must be '+' or '-', got '{cols[7]}'")
                if len(errors) >= 10:
                    errors.append("  ... (further format errors suppressed)")
                    break
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if n == 0:
        return False, ["  file is empty"]
    return len(errors) == 0, errors


def _check_sequence_lengths_format(path):
    """
    all_sequence_lengths.txt — 3 tab-separated columns:
      col0  chromosome     integer
      col1  length         integer > 0
      col2  species_label  non-empty string (e.g. 'Genus_sps1:100k')
    """
    errors = []
    n = 0
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                if lineno > FORMAT_SAMPLE_LINES:
                    break
                line = raw.rstrip("\n\r")
                if not line:
                    continue
                n += 1
                cols = line.split("\t")
                if len(cols) != 3:
                    errors.append(f"  line {lineno}: expected 3 tab-separated columns, got {len(cols)}")
                    if len(errors) >= 5:
                        errors.append("  ... (further column errors suppressed)")
                        break
                    continue
                try:
                    int(cols[0])
                except ValueError:
                    errors.append(f"  line {lineno}: col1 (chromosome) must be integer, got '{cols[0]}'")
                try:
                    ln = int(cols[1])
                    if ln <= 0:
                        errors.append(f"  line {lineno}: col2 (length={ln}) must be > 0")
                except ValueError:
                    errors.append(f"  line {lineno}: col2 (length) must be integer, got '{cols[1]}'")
                if not cols[2].strip():
                    errors.append(f"  line {lineno}: col3 (species label) is empty")
                if len(errors) >= 10:
                    errors.append("  ... (further format errors suppressed)")
                    break
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if n == 0:
        return False, ["  file is empty"]
    return len(errors) == 0, errors


def _check_fasta_format(path):
    """
    FASTA file — alternating header/sequence lines:
      • Header lines start with '>'
      • Sequence lines contain valid IUPAC nucleotide chars
    Checks first FORMAT_SAMPLE_LINES lines only for speed.
    Returns (ok, errors).
    """
    errors = []
    n_headers = 0
    n_seq_lines = 0
    in_seq = False
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                if lineno > FORMAT_SAMPLE_LINES:
                    break
                line = raw.rstrip("\n\r")
                if not line:
                    continue
                if line.startswith(">"):
                    n_headers += 1
                    in_seq = True
                    if len(line) < 2:
                        errors.append(f"  line {lineno}: header '>' has no sequence identifier")
                else:
                    if not in_seq:
                        errors.append(f"  line {lineno}: sequence line before any header '>'")
                    bad = set(line) - _VALID_NT
                    if bad:
                        errors.append(
                            f"  line {lineno}: invalid nucleotide character(s): "
                            f"{', '.join(sorted(bad))}"
                        )
                    n_seq_lines += 1
                if len(errors) >= 5:
                    errors.append("  ... (further format errors suppressed)")
                    break
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if n_headers == 0:
        return False, ["  no FASTA headers ('>') found — not a valid FASTA file"]
    if n_seq_lines == 0:
        errors.append("  no sequence lines found — FASTA appears to be header-only")
    return len(errors) == 0, errors


def _check_fasta_dir(fasta_dir, reference_name):
    """
    Checks the fasta/ directory:
      1. Reference species file present ({reference_name}.fa/.fna/.fasta)
      2. At least one query species file present (besides reference)
      3. Each FASTA file is validly formatted
    Returns (ok, errors).
    """
    errors = []
    exts = (".fa", ".fna", ".fasta")

    # Collect all FASTA files
    all_fastas = []
    for ext in exts:
        all_fastas.extend(glob.glob(os.path.join(fasta_dir, f"*{ext}")))

    if not all_fastas:
        return False, [f"  no FASTA files found in {fasta_dir}"]

    # Check reference present
    ref_found = any(
        os.path.splitext(os.path.basename(f))[0] == reference_name
        for f in all_fastas
    )
    if not ref_found:
        errors.append(
            f"  reference species file '{reference_name}.fa' (or .fna/.fasta) "
            f"not found in {fasta_dir}"
        )

    # Check at least one query species (total >= 2)
    if len(all_fastas) < 2:
        errors.append(
            f"  only {len(all_fastas)} FASTA file(s) found — need at least 1 reference + 1 query"
        )

    # Validate format of each FASTA file
    for fpath in sorted(all_fastas):
        f_ok, f_errs = _check_fasta_format(fpath)
        if not f_ok:
            errors.append(f"  {os.path.basename(fpath)} — format errors:")
            errors.extend(f_errs)

    return len(errors) == 0, errors


def _check_protein_annotation_format(path):
    """
    protein_annotation.tsv — 5 tab-separated columns:
      col0  chromosome/contig   non-empty string
      col1  start               integer >= 0
      col2  end                 integer > start
      col3  value               integer
      col4  gene/protein ID     non-empty string
    """
    errors = []
    n = 0
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                if lineno > FORMAT_SAMPLE_LINES:
                    break
                line = raw.rstrip("\n\r")
                if not line:
                    continue
                n += 1
                cols = line.split("\t")
                if len(cols) != 5:
                    errors.append(f"  line {lineno}: expected 5 tab-separated columns, got {len(cols)}")
                    if len(errors) >= 5:
                        errors.append("  ... (further column errors suppressed)")
                        break
                    continue
                # col0 chromosome — non-empty string
                if not cols[0].strip():
                    errors.append(f"  line {lineno}: col1 (chromosome) is empty")
                # col1 start — integer >= 0
                try:
                    start = int(cols[1])
                    if start < 0:
                        errors.append(f"  line {lineno}: col2 (start={start}) must be >= 0")
                except ValueError:
                    errors.append(f"  line {lineno}: col2 (start) must be integer, got '{cols[1]}'")
                # col2 end — integer > start
                try:
                    end = int(cols[2]); start_ = int(cols[1])
                    if end <= start_:
                        errors.append(f"  line {lineno}: col3 (end={end}) must be > col2 (start={start_})")
                except ValueError:
                    errors.append(f"  line {lineno}: col3 (end) must be integer, got '{cols[2]}'")
                # col3 — integer
                try:
                    int(cols[3])
                except ValueError:
                    errors.append(f"  line {lineno}: col4 must be integer, got '{cols[3]}'")
                # col4 gene ID — non-empty string
                if not cols[4].strip():
                    errors.append(f"  line {lineno}: col5 (gene/protein ID) is empty")
                if len(errors) >= 10:
                    errors.append("  ... (further format errors suppressed)")
                    break
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if n == 0:
        return False, ["  file is empty"]
    return len(errors) == 0, errors


def _check_newick_format(path):
    """
    tree.txt — Newick format:
      • Non-empty
      • Contains at least one '(' and ')'
      • Parentheses are balanced
      • Ends with ';'
    """
    errors = []
    try:
        content = open(path).read().strip()
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if not content:
        return False, ["  file is empty"]
    if not content.endswith(";"):
        errors.append("  Newick tree must end with ';'")
    if "(" not in content or ")" not in content:
        errors.append("  Newick tree must contain parentheses '(' and ')'")
    else:
        depth = 0
        for i, ch in enumerate(content):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth < 0:
                    errors.append(f"  unmatched ')' at position {i}")
                    break
        if depth != 0:
            errors.append(f"  unbalanced parentheses: {depth} '(' left unclosed")
    return len(errors) == 0, errors


def _check_species_info_format(path):
    """
    species_info.txt — 3 space-separated columns per line:
      col0  species name   non-empty string
      col1  integer        (e.g. index/order)
      col2  integer        (e.g. flag/weight)
    """
    errors = []
    n = 0
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line:
                    continue
                n += 1
                cols = line.split()
                if len(cols) != 3:
                    errors.append(
                        f"  line {lineno}: expected 3 space-separated columns, got {len(cols)}: '{line}'"
                    )
                    if len(errors) >= 5:
                        errors.append("  ... (further column errors suppressed)")
                        break
                    continue
                # col0 species name — non-empty string (already guaranteed by split())
                # col1 — integer
                try:
                    int(cols[1])
                except ValueError:
                    errors.append(f"  line {lineno}: col2 must be integer, got '{cols[1]}'")
                # col2 — integer
                try:
                    int(cols[2])
                except ValueError:
                    errors.append(f"  line {lineno}: col3 must be integer, got '{cols[2]}'")
                if len(errors) >= 10:
                    errors.append("  ... (further format errors suppressed)")
                    break
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]
    if n == 0:
        return False, ["  file is empty"]
    return len(errors) == 0, errors


def _check_axt_format(path):
    """
    AXT alignment file — records separated by blank lines.
    Each record:
      line 1 (header): 9 space-separated fields
                       alignment_id(int) chr1(str) start1(int) end1(int)
                       chr2(str) start2(int) end2(int) strand(+/-) score(int)
      line 2: reference sequence  (non-empty, IUPAC + gaps + lowercase)
      line 3: query sequence      (non-empty, IUPAC + gaps + lowercase)
    Checks first FORMAT_SAMPLE_LINES lines only.
    """
    errors = []
    valid_seq_chars = set("ACGTacgtNnRrYySsWwKkMmBbDdHhVv-")
    try:
        with open(path) as fh:
            lines = [ln.rstrip("\n\r") for ln in fh]
    except OSError as exc:
        return False, [f"  cannot read file: {exc}"]

    if not lines:
        return False, ["  file is empty"]

    # Parse records: header + seq1 + seq2 + blank
    i = 0
    record_num = 0
    checked_lines = 0
    while i < len(lines) and checked_lines < FORMAT_SAMPLE_LINES:
        # Skip blank lines between records
        if not lines[i].strip():
            i += 1
            continue

        record_num += 1
        header_lineno = i + 1

        # Header line
        header = lines[i].split()
        if len(header) != 9:
            errors.append(
                f"  record {record_num} (line {header_lineno}): AXT header must have 9 "
                f"space-separated fields, got {len(header)}: '{lines[i]}'"
            )
            i += 1
            checked_lines += 1
            if len(errors) >= 5:
                errors.append("  ... (further format errors suppressed)")
                break
            continue
        # field 0: alignment id (int)
        try:
            int(header[0])
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field1 (alignment_id) must be integer, got '{header[0]}'")
        # field 2: start1 (int >= 0)
        try:
            s1 = int(header[2])
            if s1 < 0:
                errors.append(f"  record {record_num} line {header_lineno}: field3 (start1) must be >= 0")
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field3 (start1) must be integer, got '{header[2]}'")
        # field 3: end1 (int > start1)
        try:
            e1 = int(header[3]); s1_ = int(header[2])
            if e1 <= s1_:
                errors.append(f"  record {record_num} line {header_lineno}: field4 (end1={e1}) must be > field3 (start1={s1_})")
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field4 (end1) must be integer, got '{header[3]}'")
        # field 5: start2 (int >= 0)
        try:
            s2 = int(header[5])
            if s2 < 0:
                errors.append(f"  record {record_num} line {header_lineno}: field6 (start2) must be >= 0")
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field6 (start2) must be integer, got '{header[5]}'")
        # field 6: end2 (int > start2)
        try:
            e2 = int(header[6]); s2_ = int(header[5])
            if e2 <= s2_:
                errors.append(f"  record {record_num} line {header_lineno}: field7 (end2={e2}) must be > field6 (start2={s2_})")
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field7 (end2) must be integer, got '{header[6]}'")
        # field 7: strand
        if header[7] not in ("+", "-"):
            errors.append(f"  record {record_num} line {header_lineno}: field8 (strand) must be '+' or '-', got '{header[7]}'")
        # field 8: score (int)
        try:
            int(header[8])
        except ValueError:
            errors.append(f"  record {record_num} line {header_lineno}: field9 (score) must be integer, got '{header[8]}'")

        checked_lines += 1
        i += 1

        # Sequence line 1
        if i >= len(lines) or not lines[i].strip():
            errors.append(f"  record {record_num}: missing reference sequence line")
            if len(errors) >= 5:
                errors.append("  ... (further format errors suppressed)")
                break
            continue
        seq1 = lines[i].strip()
        bad1 = set(seq1) - valid_seq_chars
        if bad1:
            errors.append(f"  record {record_num} line {i+1}: invalid char(s) in ref sequence: {', '.join(sorted(bad1))}")
        checked_lines += 1
        i += 1

        # Sequence line 2
        if i >= len(lines) or not lines[i].strip():
            errors.append(f"  record {record_num}: missing query sequence line")
            if len(errors) >= 5:
                errors.append("  ... (further format errors suppressed)")
                break
            continue
        seq2 = lines[i].strip()
        bad2 = set(seq2) - valid_seq_chars
        if bad2:
            errors.append(f"  record {record_num} line {i+1}: invalid char(s) in query sequence: {', '.join(sorted(bad2))}")
        checked_lines += 1
        i += 1

        if len(errors) >= 10:
            errors.append("  ... (further format errors suppressed)")
            break

    if record_num == 0:
        return False, ["  no AXT records found — file may be empty or malformed"]
    return len(errors) == 0, errors


def _check_axt_dir(axt_dir):
    """Validate format of every .axt file in the directory."""
    errors = []
    axt_files = sorted(glob.glob(os.path.join(axt_dir, "*.axt")))
    if not axt_files:
        return False, [f"  no .axt files found in {axt_dir}"]
    for af in axt_files:
        f_ok, f_errs = _check_axt_format(af)
        if not f_ok:
            errors.append(f"  {os.path.basename(af)}:")
            errors.extend(f_errs)
    return len(errors) == 0, errors


# ── Item class ───────────────────────────────────────────────────────────────

class Item:
    def __init__(self, label, path, kind, required_by, fmt=None, params=None):
        """
        fmt    : format-check descriptor (str or None)
        params : dict of extra parameters passed to format checker (e.g. reference_name)
        """
        self.label, self.path, self.kind = label, path, kind
        self.required_by = required_by
        self.fmt    = fmt
        self.params = params or {}

    def check(self):
        """Returns (ok: bool, detail: str, fmt_errors: list[str])"""
        p = self.path
        fmt_errors = []

        # ── Existence / type checks ──────────────────────────────────────────
        if self.kind == "file":
            ok = os.path.isfile(p)
            detail = "" if ok else f"missing file: {p}"

        elif self.kind == "dir":
            ok = os.path.isdir(p)
            detail = "" if ok else f"directory missing: {p}"

        elif self.kind == "fasta_dir":
            if not os.path.isdir(p):
                return False, f"directory missing: {p}", []
            ok, detail = True, ""

        elif self.kind == "axt_dir":
            if not os.path.isdir(p):
                return False, f"directory missing: {p}", []
            if not has_ext(p, ".axt"):
                return False, f"no .axt files found in {p}", []
            ok, detail = True, ""

        elif self.kind == "txt_dir":
            if not os.path.isdir(p):
                return False, f"directory missing: {p}", []
            if not has_ext(p, ".txt"):
                return False, f"no .txt alignment files found in {p}", []
            ok, detail = True, ""

        else:
            return False, f"unknown kind '{self.kind}'", []

        if not ok:
            return False, detail, []

        # ── Format checks (only when file/dir exists) ────────────────────────
        if self.fmt == "satsuma_alignment":
            for tf in sorted(glob.glob(os.path.join(p, "*.txt"))):
                f_ok, f_errs = _check_satsuma_format(tf)
                if not f_ok:
                    fmt_errors.append(f"  {os.path.basename(tf)}:")
                    fmt_errors.extend(f_errs)

        elif self.fmt == "sequence_lengths":
            _, fmt_errors = _check_sequence_lengths_format(p)

        elif self.fmt == "fasta_dir":
            ref = self.params.get("reference_name", "")
            _, fmt_errors = _check_fasta_dir(p, ref)

        elif self.fmt == "protein_annotation":
            _, fmt_errors = _check_protein_annotation_format(p)

        elif self.fmt == "newick":
            _, fmt_errors = _check_newick_format(p)

        elif self.fmt == "species_info":
            _, fmt_errors = _check_species_info_format(p)

        elif self.fmt == "axt_dir":
            _, fmt_errors = _check_axt_dir(p)

        return True, detail, fmt_errors


# ── Build check list ─────────────────────────────────────────────────────────

def build_checks(cfg):
    s        = cfg.get("run_stages", {})
    in_base  = cfg.get("base_input_dir", "inputs")
    ref_name = cfg.get("reference_name", "")
    items    = []

    def p(*parts):
        return os.path.join(in_base, *parts)

    run_satsuma  = s.get("run_satsuma_alignment",      False)
    run_lastz    = s.get("run_lastz_alignment",         False)
    run_synteny  = s.get("synteny_processing",          False)
    run_eba      = s.get("eba_analysis",                False)
    run_enrich   = s.get("enrichment_analysis",         False)
    run_chain    = s.get("chainNet_generation",         False)
    run_ancestor = s.get("Ancestor_seq_recunstruction", False)

    # fasta/ — existence + reference/query species + FASTA format
    if run_satsuma or run_lastz or run_chain or run_ancestor:
        fasta_dir = cfg.get("satsuma_align", {}).get("fasta_dir", p("fasta"))
        items.append(Item(
            "fasta/  (genome FASTA files)",
            fasta_dir, "fasta_dir",
            "run_satsuma_alignment / run_lastz_alignment / chainNet_generation / Ancestor_seq_recunstruction",
            fmt="fasta_dir",
            params={"reference_name": ref_name}
        ))

    # Satsuma pre-computed alignments + 8-column format
    if (not run_satsuma) and (run_synteny or run_eba or run_enrich):
        satsuma_dir = cfg.get("satsuma_alignments",
                              p("synteny_processing", "Satsuma_alignments"))
        items.append(Item(
            "synteny_processing/Satsuma_alignments/  (pre-computed .txt alignment files)",
            satsuma_dir, "txt_dir",
            "synteny_processing / eba_analysis / enrichment_analysis  (run_satsuma_alignment: false)",
            fmt="satsuma_alignment"
        ))

    # all_sequence_lengths.txt + 3-column format
    # This file is optional when FASTA files are present — the pipeline will
    # auto-generate it via genome_length_maker.sh before synteny processing.
    if run_synteny:
        seq_len = cfg.get("sequence_lengths_file",
                          p("synteny_processing", "all_sequence_lengths.txt"))
        fasta_dir = cfg.get("satsuma_align", {}).get("fasta_dir", p("fasta"))
        _has_fastas = os.path.isdir(fasta_dir) and any(
            glob.glob(os.path.join(fasta_dir, f"*{ext}"))
            for ext in (".fa", ".fna", ".fasta")
        )
        if os.path.isfile(seq_len):
            # File already exists — validate its format
            items.append(Item(
                "synteny_processing/all_sequence_lengths.txt",
                seq_len, "file", "synteny_processing",
                fmt="sequence_lengths"
            ))
        elif _has_fastas:
            # File missing but FASTAs are present → will be auto-generated; skip hard check
            print(f"  {CYAN}ℹ{NC}  synteny_processing/all_sequence_lengths.txt not found "
                  f"— will be auto-generated from FASTA files via genome_length_maker.sh")
        else:
            # Neither the file nor FASTAs exist → hard failure
            items.append(Item(
                "synteny_processing/all_sequence_lengths.txt",
                seq_len, "file", "synteny_processing",
                fmt="sequence_lengths"
            ))

    # classification.eba — existence only (EBA-tool-specific binary/text format)
    if run_eba:
        eba_c = cfg.get("eba", {}).get("c", p("eba_analysis", "classification.eba"))
        items.append(Item(
            "eba_analysis/classification.eba",
            eba_c, "file", "eba_analysis"
        ))

    # protein_annotation.tsv + 5-column format
    if run_enrich:
        ann = cfg.get("enrichment", {}).get("annotation_file",
              p("enrichment_analysis", "protein_annotation.tsv"))
        items.append(Item(
            "enrichment_analysis/protein_annotation.tsv",
            ann, "file", "enrichment_analysis",
            fmt="protein_annotation"
        ))

    # Pre-computed .axt files + AXT format
    if (not run_lastz) and (run_chain or run_ancestor):
        axt_dir = cfg.get("chainNet", {}).get("lastZ_alignments",
                  p("Ancestor_seq_recunstruction", "LastZ_alignments"))
        items.append(Item(
            "Ancestor_seq_recunstruction/LastZ_alignments/  (pre-computed .axt files)",
            axt_dir, "axt_dir",
            "chainNet_generation / Ancestor_seq_recunstruction  (run_lastz_alignment: false)",
            fmt="axt_dir"
        ))

    # tree.txt + Newick format
    if run_ancestor:
        tree = cfg.get("deschrambler", {}).get("tree_file",
               p("Ancestor_seq_recunstruction", "tree.txt"))
        items.append(Item(
            "Ancestor_seq_recunstruction/tree.txt",
            tree, "file", "Ancestor_seq_recunstruction",
            fmt="newick"
        ))

    # species_info.txt + 3-column space-separated format
    if run_ancestor or run_chain:
        spp = cfg.get("deschrambler", {}).get("species_file",
              p("Ancestor_seq_recunstruction", "species_info.txt"))
        items.append(Item(
            "Ancestor_seq_recunstruction/species_info.txt",
            spp, "file", "Ancestor_seq_recunstruction / chainNet_generation",
            fmt="species_info"
        ))

    return items


# ── Report printer ───────────────────────────────────────────────────────────

def print_report(cfg, items, verbose):
    s       = cfg.get("run_stages", {})
    in_base = cfg.get("base_input_dir", "inputs")
    active  = [k for k, v in s.items() if v]

    print()
    print(f"{BOLD}╔══════════════════════════════════════════════════════════╗{NC}")
    print(f"{BOLD}║           SYBR Pipeline — Input Verification             ║{NC}")
    print(f"{BOLD}╚══════════════════════════════════════════════════════════╝{NC}")
    print(f"Active stages : {CYAN}{', '.join(active) if active else '(none)'}{NC}")
    print(f"Input base dir: {in_base}")
    print(f"┌─ Existence & format checks ──────────────────────────────┐")

    exist_errors  = 0
    format_errors = 0

    for item in items:
        ok, detail, fmt_errs = item.check()
        if not ok:
            exist_errors += 1
            print(f"  |-- {CROSS} {item.label}")
            print(f"      {RED}[{detail}]{NC}")
        elif fmt_errs:
            format_errors += 1
            print(f"  |-- {WARN} {item.label}  {YELLOW}[exists, format errors below]{NC}")
            for fe in fmt_errs:
                print(f"      {RED}{fe}{NC}")
        else:
            fmt_note = f"  {GREEN}(format OK){NC}" if item.fmt else ""
            print(f"  |-- {CHECK} {item.label}{fmt_note}")

    if not items:
        print(f"  (no input checks required for active modules)")

    total_errors = exist_errors + format_errors
    result_str = (
        "PASSED" if total_errors == 0
        else f"FAILED ({exist_errors} missing, {format_errors} format error(s))"
    )
    print(f"└─ Result: {result_str} {'─' * 10}┘")
    print()

    if total_errors > 0:
        if exist_errors:
            print(f"{RED}  {exist_errors} required file(s)/directory(ies) are missing.{NC}")
        if format_errors:
            print(f"{RED}  {format_errors} file(s) have format errors — see details above.{NC}")
        print(f"Resolve the issues above before re-running.")
        print(f"{'═' * 60}")
        print(f"{RED}✗ Pre-flight check failed — {total_errors} issue(s) found.{NC}")
        print(f"{'═' * 60}")
    else:
        print(f"{'═' * 60}")
        print(f"{GREEN}✔ All required inputs present and correctly formatted.{NC}")
        print(f"{'═' * 60}")

    return total_errors


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="SYBR pipeline input validator")
    parser.add_argument("--config",  required=True, help="Merged YAML config file")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        print(f"ERROR: config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    cfg    = load_config(args.config)
    items  = build_checks(cfg)
    errors = print_report(cfg, items, args.verbose)
    sys.exit(0 if errors == 0 else 1)


if __name__ == "__main__":
    main()
