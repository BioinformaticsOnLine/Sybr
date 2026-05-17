#!/usr/bin/env Rscript
# linear_synteny.R
#
# Usage: Rscript linear_synteny.R <lengths_file> <split_by_species_dir>
#
# The full all_sequence_lengths.txt is passed in. This script filters it to
# only the entries relevant to the current resolution (derived from the
# split_by_species_dir path, e.g. ".../100k/split_by_species" → "100k").
#
# Column conventions:
#   Lengths file : chr_id, length, species_label   (col3 may have ":100k" suffix)
#   Pipeline aln : col8=reference_label, col9=query_label
#   draw.linear  : col8=query_label,     col9=reference_label  → swap cols 8 & 9

suppressPackageStartupMessages(library(syntenyPlotteR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript linear_synteny.R <lengths_file> <split_by_species_dir>")

lengths_file <- normalizePath(args[1])
split_dir    <- normalizePath(args[2])

if (!file.exists(lengths_file)) stop(paste("Lengths file not found:", lengths_file))
if (!dir.exists(split_dir))     stop(paste("Directory not found:", split_dir))

pdf(NULL)

# ── Derive resolution tag from directory path ─────────────────────────────────
# e.g.  .../synteny_plots/100k/split_by_species  →  res = "100k"
res_name <- basename(dirname(split_dir))
cat("Resolution   :", res_name, "\n")
cat("Lengths file :", lengths_file, "\n")
cat("Directory    :", split_dir, "\n")

# ── Load and filter lengths file ──────────────────────────────────────────────
# Keep only rows for this resolution (col3 ends with ":100k") or query species
# (col3 has no ":" — no resolution suffix).
ldat <- read.table(lengths_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

# Filter: reference rows for this resolution OR query rows (no ":")
keep <- grepl(paste0(":", res_name, "$"), ldat[[3]]) | !grepl(":", ldat[[3]])
ldat <- ldat[keep, , drop=FALSE]

# Strip resolution suffix from reference label so names match alignment files
ldat[[3]] <- sub(":.*$", "", ldat[[3]])

cat("Lengths rows after filtering:", nrow(ldat), "\n")
cat("Species in lengths:\n")
for (s in unique(ldat[[3]])) cat(" -", s, "\n")

# ── Collect alignment files ───────────────────────────────────────────────────
alignment_files <- list.files(split_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(alignment_files) == 0) stop(paste("No .txt files found in", split_dir))
cat("Alignment files:\n")
for (f in alignment_files) cat(" -", basename(f), "\n")

# ── Write temp lengths file ───────────────────────────────────────────────────
tmp_lengths <- tempfile(fileext = ".txt")
write.table(ldat, tmp_lengths, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
on.exit(unlink(tmp_lengths), add=TRUE)

# ── Build cleaned temp alignment files ───────────────────────────────────────
# Swap col8/col9: pipeline has (ref, query), draw.linear wants (query, ref)
# Strip resolution suffix from col9 (reference) after swap
clean_align_files <- vapply(alignment_files, function(f) {
  dat      <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  tmp8     <- dat[[8]]
  dat[[8]] <- dat[[9]]                 # query  → col8
  dat[[9]] <- sub(":.*$", "", tmp8)    # ref    → col9, suffix stripped
  tmp <- tempfile(fileext=".txt")
  write.table(dat, tmp, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  tmp
}, character(1))
on.exit(unlink(clean_align_files), add=TRUE)

# ── Call draw.linear once per alignment file ──────────────────────────────────
for (i in seq_along(alignment_files)) {
  input_file  <- alignment_files[[i]]
  clean_file  <- clean_align_files[[i]]
  output_name <- tools::file_path_sans_ext(basename(input_file))

  cat("\nProcessing:", basename(input_file), "→", paste0(output_name, ".pdf"), "\n")

  tryCatch({
    draw.linear(output_name, tmp_lengths, clean_file,
                fileformat = "pdf", w = 13, h = 5, opacity = 0.5,
                directory = split_dir)
    cat("Done:", file.path(split_dir, paste0(output_name, ".pdf")), "\n")
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}

while (!is.null(dev.list())) dev.off()
cat("Done.\n")
