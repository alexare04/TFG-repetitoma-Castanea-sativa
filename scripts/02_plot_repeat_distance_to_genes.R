# ============================================================
# Script: 10_plot_repeat_distance_to_genes.R
# Purpose: Plot the distance of non-overlapping repetitive
#          elements to the nearest annotated gene, grouped by
#          repeat class.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(tidyverse)

# -----------------------------
# Input files
# -----------------------------

distance_file <- "results/gene_overlap_distances/RM_geneImpact_grouped/distances_by_group.long.tsv"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

summary_output <- "results/tables/repeat_distance_to_genes_summary.tsv"
plot_output <- "results/figures/repeat_distance_to_genes_by_class.pdf"

# -----------------------------
# Load data
# -----------------------------

dist_long <- read_tsv(
  distance_file,
  show_col_types = FALSE
)

# -----------------------------
# Filter non-overlapping repeats
# -----------------------------

# Distance = 0 corresponds to repetitive elements that overlap genes.
# For this plot, only non-overlapping repeats are represented.

dist_long_nonzero <- dist_long %>%
  filter(distance_bp > 0)

# -----------------------------
# Calculate summary statistics
# -----------------------------

dist_sum_nonzero <- dist_long_nonzero %>%
  group_by(group) %>%
  summarise(
    n = n(),
    min = min(distance_bp, na.rm = TRUE),
    q1 = quantile(distance_bp, 0.25, na.rm = TRUE),
    median = median(distance_bp, na.rm = TRUE),
    mean = mean(distance_bp, na.rm = TRUE),
    q3 = quantile(distance_bp, 0.75, na.rm = TRUE),
    max = max(distance_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median)

write.table(
  dist_sum_nonzero,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(dist_sum_nonzero)

# -----------------------------
# Order repeat classes by median distance
# -----------------------------

dist_long_nonzero <- dist_long_nonzero %>%
  mutate(
    group = factor(group, levels = dist_sum_nonzero$group)
  )

# -----------------------------
# Plot distance to nearest gene
# -----------------------------

p <- ggplot(dist_long_nonzero, aes(x = group, y = distance_bp)) +
  geom_violin(fill = "grey80", color = "grey30", trim = TRUE) +
  coord_flip() +
  scale_y_continuous(
    trans = "log10",
    breaks = c(10, 100, 1000, 10000, 100000, 300000),
    labels = c("10", "100", "1 kb", "10 kb", "100 kb", "300 kb")
  ) +
  labs(
    x = NULL,
    y = "Distance to nearest gene (bp, non-overlapping repeats only)",
    title = "Distance of non-overlapping repetitive elements to nearest gene"
  ) +
  theme_minimal()

pdf(plot_output, width = 10, height = 6)
print(p)
dev.off()
