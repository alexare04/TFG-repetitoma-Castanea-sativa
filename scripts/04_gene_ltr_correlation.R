# ============================================================
# Script: 04_gene_ltr_correlation.R
# Purpose: Calculate Spearman correlations between gene density
#          and Gypsy/Copia density using 1 Mb windows.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# Input files
# -----------------------------

window_counts_file <- "results/tables/window_counts_1Mb.tsv"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

global_output <- "results/tables/global_spearman_gene_ltr.tsv"
chrom_output <- "results/tables/chromosome_spearman_gene_ltr.tsv"
plot_output <- "results/figures/Gypsy_Copia_gene_density_correlations.pdf"

# -----------------------------
# Load data
# -----------------------------

df_all <- read_tsv(window_counts_file, show_col_types = FALSE)

# Remove windows without gene annotations
df_all <- df_all %>%
  filter(genes > 0)

# -----------------------------
# Global Spearman correlations
# -----------------------------

cor_gypsy <- cor.test(df_all$genes, df_all$gypsy, method = "spearman")
cor_copia <- cor.test(df_all$genes, df_all$copia, method = "spearman")

global_cor <- data.frame(
  comparison = c("Gypsy_vs_genes", "Copia_vs_genes"),
  rho = c(unname(cor_gypsy$estimate), unname(cor_copia$estimate)),
  p_value = c(cor_gypsy$p.value, cor_copia$p.value)
)

write.table(
  global_cor,
  file = global_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Spearman correlations by chromosome
# -----------------------------

cor_chr <- df_all %>%
  group_by(seqnames) %>%
  summarise(
    n_windows = n(),
    rho_gypsy = cor(genes, gypsy, method = "spearman"),
    p_gypsy = cor.test(genes, gypsy, method = "spearman")$p.value,
    rho_copia = cor(genes, copia, method = "spearman"),
    p_copia = cor.test(genes, copia, method = "spearman")$p.value,
    .groups = "drop"
  )

write.table(
  cor_chr,
  file = chrom_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Scatter plots
# -----------------------------

p_gypsy <- ggplot(df_all, aes(x = genes, y = gypsy)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_smooth(method = "lm", color = "black") +
  coord_cartesian(ylim = c(0, 200)) +
  labs(
    x = "Genes / 1 Mb",
    y = "LTR retrotransposons / 1 Mb",
    title = "Relationship between Gypsy density and gene density"
  ) +
  theme_minimal()

p_copia <- ggplot(df_all, aes(x = genes, y = copia)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_smooth(method = "lm", color = "black") +
  coord_cartesian(ylim = c(0, 200)) +
  labs(
    x = "Genes / 1 Mb",
    y = "LTR retrotransposons / 1 Mb",
    title = "Relationship between Copia density and gene density"
  ) +
  theme_minimal()

pdf(plot_output, width = 10, height = 5)
print(p_gypsy)
print(p_copia)
dev.off()
