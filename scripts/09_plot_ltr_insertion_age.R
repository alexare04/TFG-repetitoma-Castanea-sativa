# ============================================================
# Script: 09_plot_ltr_insertion_age.R
# Purpose: Plot insertion age distributions of intact LTR
#          retrotransposons located inside and outside genes,
#          and calculate summary statistics by family and location.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(ggplot2)
library(dplyr)
library(readr)

# -----------------------------
# Input files
# -----------------------------

gypsy_in_file <- "results/ltr_insertion_age/Gypsy_in_genes_age_MYA.txt"
copia_in_file <- "results/ltr_insertion_age/Copia_in_genes_age_MYA.txt"
gypsy_out_file <- "results/ltr_insertion_age/Gypsy_outside_genes_age_MYA.txt"
copia_out_file <- "results/ltr_insertion_age/Copia_outside_genes_age_MYA.txt"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

summary_output <- "results/tables/ltr_insertion_age_summary.tsv"
plot_output <- "results/figures/LTR_insertion_age_inside_outside_genes.pdf"

# -----------------------------
# Load data
# -----------------------------

gypsy_in <- read_tsv(
  gypsy_in_file,
  col_names = "age",
  show_col_types = FALSE
) %>%
  mutate(
    family = "Gypsy",
    location = "Inside genes"
  )

copia_in <- read_tsv(
  copia_in_file,
  col_names = "age",
  show_col_types = FALSE
) %>%
  mutate(
    family = "Copia",
    location = "Inside genes"
  )

gypsy_out <- read_tsv(
  gypsy_out_file,
  col_names = "age",
  show_col_types = FALSE
) %>%
  mutate(
    family = "Gypsy",
    location = "Outside genes"
  )

copia_out <- read_tsv(
  copia_out_file,
  col_names = "age",
  show_col_types = FALSE
) %>%
  mutate(
    family = "Copia",
    location = "Outside genes"
  )

df_all <- bind_rows(
  gypsy_in,
  copia_in,
  gypsy_out,
  copia_out
)

df_all$family <- factor(df_all$family, levels = c("Copia", "Gypsy"))
df_all$location <- factor(df_all$location, levels = c("Inside genes", "Outside genes"))

# -----------------------------
# Summary statistics
# -----------------------------

age_summary <- df_all %>%
  group_by(family, location) %>%
  summarise(
    n = n(),
    median_age = median(age, na.rm = TRUE),
    mean_age = mean(age, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  age_summary,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(age_summary)

# -----------------------------
# Plot insertion age distributions
# -----------------------------

p <- ggplot(df_all, aes(x = age, color = family)) +
  geom_density(linewidth = 1.2) +
  facet_wrap(~location) +
  theme_minimal() +
  labs(
    x = "Estimated age (million years)",
    y = "Density",
    color = "LTR Superfamily"
  )

pdf(plot_output, width = 10, height = 5)
print(p)
dev.off()
