#!/usr/bin/env Rscript

# Load library
suppressPackageStartupMessages(library(syntenyPlotteR))

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript run_synteny_linear.R <lengths_file> <parent_folder>")
}

lengths_file   <- args[1]
parent_folder  <- args[2]

# Prevent creation of Rplots.pdf by using a null device
pdf(NULL)

# Function to find all alignment files in split_by_species subdirectories
find_alignment_files <- function(parent_folder) {
  # Get all directories under parent_folder (should be 100, 300, 500, etc.)
  subdirs <- list.dirs(parent_folder, recursive = FALSE)
  
  # Find split_by_species directories within each subdirectory
  alignment_files <- list()
  
  for (subdir in subdirs) {
    split_dir <- file.path(subdir, "split_by_species")
    if (dir.exists(split_dir)) {
      # Get all .txt files in the split_by_species directory
      files <- list.files(split_dir, pattern = "\\.txt$", full.names = TRUE)
      alignment_files <- c(alignment_files, files)
    }
  }
  
  return(alignment_files)
}

# Find all alignment files
alignment_files <- find_alignment_files(parent_folder)

if (length(alignment_files) == 0) {
  stop(paste("No alignment files found in split_by_species directories under", parent_folder))
}

cat("Found", length(alignment_files), "alignment files to process:\n")
for (file in alignment_files) {
  cat("-", file, "\n")
}

# Process each alignment file
for (alignment_file in alignment_files) {
  # Determine directory of input file
  directory <- dirname(normalizePath(alignment_file))
  
  # Derive output name from alignment file (remove extension)
  output_name <- tools::file_path_sans_ext(basename(alignment_file))
  
  cat("\nProcessing:", alignment_file, "\n")
  cat("Output will be saved to:", directory, "\n")
  cat("Output name:", output_name, "\n")
  
  # Run syntenyPlotteR
  tryCatch({
    draw.linear(
      output_name,
      lengths_file,
      alignment_file,
      directory = directory
    )
    cat("Successfully processed:", alignment_file, "\n")
  }, error = function(e) {
    cat("Error processing", alignment_file, ":", e$message, "\n")
  })
}

# Ensure any remaining graphics devices are closed
while (!is.null(dev.list())) {
  dev.off()
}

cat("\nAll files processed!\n")
