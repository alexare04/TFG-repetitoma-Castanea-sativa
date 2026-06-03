# ============================================================
# Script: 03_ltr_distribution_windows.R
# Purpose: Calculate gene, Gypsy and Copia density in 1 Mb
#          genomic windows and generate chromosome-wise plots.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(Biostrings)
library(ggplot2)
library(patchwork)
library(dplyr)

# -----------------------------
# Input files
# -----------------------------

genes_file <- "data/genomic.gff"
ltr_file   <- "data/castano_genoma_filtrado.fasta.out.gff3"
genome_file <- "data/castano_genoma_filtrado.fasta"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

output_pdf <- "results/figures/LTR_Gypsy_Copia_Genes_by_chromosome.pdf"
output_table <- "results/tables/window_counts_1Mb.tsv"

# -----------------------------
# Import data
# -----------------------------

genes <- import(genes_file)
ltr <- import(ltr_file)
genome <- readDNAStringSet(genome_file)

# Keep only gene features
genes <- genes[genes$type == "gene"]

# Filter Gypsy and Copia elements
gypsy <- ltr[grepl("LTR/Gypsy", ltr$Classification)]
copia <- ltr[grepl("LTR/Copia", ltr$Classification)]

# -----------------------------
# Harmonize chromosome names
# -----------------------------

common_seqlevels <- intersect(seqlevels(genes), names(genome))

genes <- keepSeqlevels(genes, common_seqlevels, pruning.mode = "coarse")
gypsy <- keepSeqlevels(gypsy, common_seqlevels, pruning.mode = "coarse")
copia <- keepSeqlevels(copia, common_seqlevels, pruning.mode = "coarse")

chrom_lengths <- setNames(width(genome[common_seqlevels]), common_seqlevels)

seqlengths(genes) <- chrom_lengths
seqlengths(gypsy) <- chrom_lengths
seqlengths(copia) <- chrom_lengths

# -----------------------------
# Create 1 Mb windows
# -----------------------------

windows <- tileGenome(
  seqlengths(genes),
  tilewidth = 1e6,
  cut.last.tile.in.chrom = TRUE
)

# Count overlaps per window
windows$genes <- countOverlaps(windows, genes)
windows$gypsy <- countOverlaps(windows, gypsy)
windows$copia <- countOverlaps(windows, copia)

# Save window counts
windows_df <- as.data.frame(windows)
write.table(
  windows_df,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Plot chromosome-wise density
# -----------------------------

chroms <- seqlevels(genes)

pdf(output_pdf, width = 10, height = 8)

for (chr in chroms) {
  
  df_chr <- as.data.frame(windows[seqnames(windows) == chr])
  df_chr$mid <- (df_chr$start + df_chr$end) / 2 / 1e6
  
  p_gypsy <- ggplot(df_chr, aes(x = mid, y = gypsy)) +
    geom_area(alpha = 0.8) +
    labs(y = "Gypsy / 1 Mb", title = chr) +
    theme_minimal()
  
  p_copia <- ggplot(df_chr, aes(x = mid, y = copia)) +
    geom_area(alpha = 0.8) +
    labs(y = "Copia / 1 Mb", x = NULL) +
    theme_minimal()
  
  p_genes <- ggplot(df_chr, aes(x = mid)) +
    geom_tile(
      aes(y = 1, fill = genes),
      width = 1,
      height = 1
    ) +
    scale_fill_gradient(low = "white", high = "black") +
    labs(
      x = "Position (Mb)",
      y = "Genes / 1 Mb",
      fill = "Genes"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  print(p_gypsy / p_copia / p_genes)
}

dev.off()
