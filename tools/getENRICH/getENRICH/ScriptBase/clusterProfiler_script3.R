library(jsonlite)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(ontologyIndex)

# Get the present working directory
pwd <- getwd()

# Load configuration
config <- fromJSON("config.json")

# Extract paths from configuration and convert to absolute paths
outdir <- file.path(config$output_files$outdir)
graph <- file.path(outdir, config$output_files$graph)
setwd(graph)

# Shared output file
enrichment_KEGG_results_csv <- file.path(outdir, config$output_files$enrichment_KEGG_results_csv)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
pvalue_threshold <- as.numeric(args[1])
organism <- args[2]

if (organism == "ko") {
  kegg_annotationTOgenes_sb_3 <- file.path(config$input_files$kegg_annotationTOgenes_sb_3)
  background_genes_sb <- file.path(config$input_files$background_genes_sb_1)
  genes_of_interest_sb <- file.path(config$input_files$genes_of_interest_sb_2)

  eggNOG_kegg <- read_tsv(kegg_annotationTOgenes_sb_3)
  background_genes <- read_tsv(background_genes_sb) %>% unlist() %>% as.vector()
  interesting_set <- read_tsv(genes_of_interest_sb) %>% unlist() %>% as.vector()

  background_kegg <- eggNOG_kegg %>%
    filter(gene %in% background_genes) %>%
    select(term) %>%
    unlist() %>%
    as.vector()

  interesting_set_kegg <- eggNOG_kegg %>%
    filter(gene %in% interesting_set) %>%
    select(term) %>%
    unlist() %>%
    as.vector()

  enrichment_kegg <- enrichKEGG(
    interesting_set_kegg,
    organism = organism,
    keyType = "kegg",
    pvalueCutoff = pvalue_threshold,
    pAdjustMethod = "BH",
    universe = background_kegg,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.05,
    use_internal_data = FALSE
  )

} else if (organism == "go") {
  # Load GO-specific inputs
  term2gene_path <- file.path(config$input_files$go_term2gene)
  go_obo_path <- file.path(config$input_files$go_obo)
  background_genes_path <- file.path(config$input_files$background_genes_sb_1)
  interesting_genes_path <- file.path(config$input_files$genes_of_interest_sb_2)

  term2gene <- read_tsv(term2gene_path, show_col_types = FALSE)

  ontology <- get_ontology(
    file = go_obo_path,
    propagate_relationships = "is_a",
    extract_tags = "everything",
    merge_equivalent_terms = TRUE
  )

  term2name <- term2gene %>%
    mutate(name = ontology$name[term]) %>%
    select(term, name) %>%
    distinct() %>%
    drop_na() %>%
    filter(!grepl("obsolete", name))

  term2gene <- term2gene %>%
    filter(term %in% term2name$term)

  background_genes <- read_tsv(background_genes_path, show_col_types = FALSE) %>%
    select(geneID) %>%
    unlist() %>%
    as.vector()

  interesting_set <- read_tsv(interesting_genes_path, show_col_types = FALSE) %>%
    filter(abs(logFC) >= 1 & adj.P.Val <= 0.05) %>%
    select(geneID) %>%
    unlist() %>%
    as.vector()

  enrichment_kegg <- enricher(
    interesting_set,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = pvalue_threshold,
    universe = background_genes,
    qvalueCutoff = 0.05
  )

} else {
  # Assume organism is a supported KEGG organism (e.g., "hsa")
  background_genes_sb <- file.path(config$input_files$background_genes_sb_1)
  genes_of_interest_sb <- file.path(config$input_files$genes_of_interest_sb_2)

  background <- read_tsv(background_genes_sb) %>% unlist() %>% as.vector()
  foreground <- read_tsv(genes_of_interest_sb) %>% unlist() %>% as.vector()

  enrichment_kegg <- enrichKEGG(
    gene = foreground,
    organism = organism,
    keyType = "ncbi-geneid",
    pvalueCutoff = pvalue_threshold,
    pAdjustMethod = "BH",
    universe = background,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.05
  )
}

# Save the result
write.csv(file = enrichment_KEGG_results_csv,
          x = enrichment_kegg@result)

message("Enrichment analysis for ", organism, " completed and saved to: ", enrichment_KEGG_results_csv)
