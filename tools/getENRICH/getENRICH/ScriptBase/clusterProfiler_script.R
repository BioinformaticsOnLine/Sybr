library(jsonlite)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(ontologyIndex)

# Get the present working directory
pwd <- getwd()

# Load configuration
#config <- fromJSON("config.json")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]
pvalue_threshold <- as.numeric(args[2])
organism <- args[3]

# Load configuration from provided path
config <- fromJSON(config_path)


# Extract paths from configuration and convert to absolute paths
outdir <- file.path(config$output_files$outdir)
graph <- file.path(outdir, config$output_files$graph)
setwd(graph)

# Shared output file
enrichment_KEGG_results_csv <- file.path(outdir, config$output_files$enrichment_KEGG_results_csv)


# Check if the organism is the default "ko" or a specified one
if (organism == "ko") {
  # Load KO-specific inputs
  kegg_annotationTOgenes_sb_3 <- file.path(config$input_files$kegg_annotationTOgenes_sb_3)
  background_genes_sb <- file.path(config$input_files$background_genes_sb_1)
  genes_of_interest_sb <- file.path(config$input_files$genes_of_interest_sb_2)
  term2gene <- file.path(config$db$term2gene)
  term2name <- file.path(config$db$term2name)

  # KEGG enrichment with KEGG orthologs
  eggNOG_kegg <- read_tsv(kegg_annotationTOgenes_sb_3)

  background_genes <- read_tsv(background_genes_sb) %>%
    unlist() %>%
    as.vector()

  # read the gene list of interest
  interesting_set <- read_tsv(genes_of_interest_sb) %>%
    unlist() %>%
    as.vector()

  # Create a clean list of KEGG orthologs for the background
  background_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% background_genes) %>%
    dplyr::select(term) %>%  # Select only the KEGG ortholog column
    unlist() %>%
    as.vector()

  # Create a clean list of KEGG orthologs for the genes of interest
  interesting_set_kegg <- eggNOG_kegg %>%
    dplyr::filter(gene %in% interesting_set) %>%
    dplyr::select(term) %>%  # Select only the KEGG ortholog column
    unlist() %>%
    as.vector()


# Read term2gene and term2name files
  term2gene_df <- read.delim(term2gene, header = FALSE)
  term2name_df <- read.delim(term2name, header = FALSE)
  
  # Format for clusterProfiler
  term2gene_formatted <- data.frame(term = term2gene_df[,1], gene = term2gene_df[,2])
  term2name_formatted <- data.frame(term = term2name_df[,1], name = term2name_df[,2])


  # Use enricher instead of enrichKEGG to allow custom term mappings
  enrichment_kegg <- enricher(gene = interesting_set_kegg,
                        TERM2GENE = term2gene_formatted,
                        TERM2NAME = term2name_formatted,
                        pvalueCutoff = pvalue_threshold,
                        pAdjustMethod = "BH",
                        universe = background_kegg,
                        minGSSize = 10,
                        maxGSSize = 500,
                        qvalueCutoff = 0.05) 



} else if (organism == "go") {
  # Load GO-specific inputs
  term2gene_path <- file.path(config$input_files$go_term2gene)
  go_obo_path <- file.path(config$db$go_obo)
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

  
  # Read preprocessed background and interesting gene lists
  background_genes <- readLines(background_genes_path)
  interesting_set <- readLines(interesting_genes_path)

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
