# ============================================================
# Script: 11_go_enrichment_ltr_genes.R
# Purpose: Perform GO enrichment analysis for genes associated
#          with Gypsy and Copia LTR retrotransposons.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(clusterProfiler)
library(readr)
library(dplyr)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)
library(patchwork)

# -----------------------------
# Input files
# -----------------------------

genes_all_file <- "data/all_genes_mappable.txt"
genes_gypsy_file <- "data/genes_with_Gypsy_mappable.txt"
genes_copia_file <- "data/genes_with_Copia_mappable.txt"
gene2go_file <- "data/gene2go.tsv"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

gypsy_table_out <- "results/tables/GO_enrichment_Gypsy.tsv"
copia_table_out <- "results/tables/GO_enrichment_Copia.tsv"
plot_out <- "results/figures/GO_enrichment_Gypsy_Copia.pdf"

# -----------------------------
# Load data
# -----------------------------

genes_all <- read_lines(genes_all_file)
genes_gypsy <- read_lines(genes_gypsy_file)
genes_copia <- read_lines(genes_copia_file)

gene2go <- read.delim(
  gene2go_file,
  header = FALSE,
  sep = "\t",
  col.names = c("gene", "GO"),
  stringsAsFactors = FALSE
)

term2gene <- gene2go[, c("GO", "gene")]

# -----------------------------
# GO enrichment
# -----------------------------

ego_gypsy <- enricher(
  gene = genes_gypsy,
  universe = genes_all,
  TERM2GENE = term2gene
)

ego_copia <- enricher(
  gene = genes_copia,
  universe = genes_all,
  TERM2GENE = term2gene
)

gypsy_df <- as.data.frame(ego_gypsy)
copia_df <- as.data.frame(ego_copia)

# -----------------------------
# Add GO term names
# -----------------------------

go_names_gypsy <- AnnotationDbi::select(
  GO.db,
  keys = gypsy_df$ID,
  columns = "TERM",
  keytype = "GOID"
)

go_names_copia <- AnnotationDbi::select(
  GO.db,
  keys = copia_df$ID,
  columns = "TERM",
  keytype = "GOID"
)

gypsy_df <- merge(
  gypsy_df,
  go_names_gypsy,
  by.x = "ID",
  by.y = "GOID",
  all.x = TRUE
)

copia_df <- merge(
  copia_df,
  go_names_copia,
  by.x = "ID",
  by.y = "GOID",
  all.x = TRUE
)

gypsy_df$Description <- gypsy_df$TERM
copia_df$Description <- copia_df$TERM

gypsy_df <- gypsy_df[!is.na(gypsy_df$Description), ]
copia_df <- copia_df[!is.na(copia_df$Description), ]

# -----------------------------
# Save tables
# -----------------------------

write.table(
  gypsy_df,
  file = gypsy_table_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  copia_df,
  file = copia_table_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Prepare plots
# -----------------------------

gypsy_plot <- gypsy_df %>%
  arrange(p.adjust) %>%
  head(20)

copia_plot <- copia_df %>%
  arrange(p.adjust) %>%
  head(20)

p_gypsy <- ggplot(
  gypsy_plot,
  aes(
    x = FoldEnrichment,
    y = reorder(Description, FoldEnrichment),
    size = Count,
    color = p.adjust
  )
) +
  geom_point() +
  xlim(1, 6) +
  theme_minimal() +
  labs(
    title = "Gypsy",
    x = "Fold enrichment",
    y = "GO term",
    color = "Adjusted p-value",
    size = "Count"
  )

p_copia <- ggplot(
  copia_plot,
  aes(
    x = FoldEnrichment,
    y = reorder(Description, FoldEnrichment),
    size = Count,
    color = p.adjust
  )
) +
  geom_point() +
  xlim(1, 6) +
  theme_minimal() +
  labs(
    title = "Copia",
    x = "Fold enrichment",
    y = "GO term",
    color = "Adjusted p-value",
    size = "Count"
  )

pdf(plot_out, width = 10, height = 12)
print(p_gypsy / p_copia)
dev.off()
