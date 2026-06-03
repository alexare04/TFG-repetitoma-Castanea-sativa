# ============================================================
# Script: 12_structural_variants_repeat_density.R
# Purpose: Analyze repetitive element density inside and outside
#          large structural variant regions between haplotypes.
# Author: Alejandro Arévalo Sánchez
# ============================================================

# -----------------------------
# Load libraries
# -----------------------------

library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(dplyr)
library(ggplot2)

# -----------------------------
# Input files
# -----------------------------

repeatmasker_gff <- "data/castano_genoma_filtrado.fasta.mod.out.gff"

# -----------------------------
# Output files
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

summary_output <- "results/tables/SV_repeat_density_summary.tsv"
wilcox_output <- "results/tables/SV_repeat_density_wilcoxon.tsv"
plot_output <- "results/figures/SV_repeat_density_by_chromosome.pdf"

# -----------------------------
# Define structural variant regions
# -----------------------------

sv_regions <- data.frame(
  chr = c("NC_134017.1", "NC_134019.1", "NC_134020.1"),
  start = c(3644689, 11157318, 1),
  end = c(6474376, 11672496, 392251),
  sv_id = c("chr5_SV", "chr7_SV", "chr8_SV")
)

# Chromosome lengths for chromosomes with SVs
chrom_lengths <- c(
  "NC_134017.1" = 81892547,
  "NC_134019.1" = 47783066,
  "NC_134020.1" = 56333209
)

# -----------------------------
# Create genomic windows
# -----------------------------

windows <- tileGenome(
  chrom_lengths,
  tilewidth = 5e5,
  cut.last.tile.in.chrom = TRUE
)

# -----------------------------
# Import repetitive elements
# -----------------------------

rep <- import(repeatmasker_gff)
rep <- keepSeqlevels(rep, names(chrom_lengths), pruning.mode = "coarse")
seqlengths(rep) <- chrom_lengths

# -----------------------------
# Calculate repeat coverage per window
# -----------------------------

hits <- findOverlaps(windows, rep)

ov <- pintersect(
  windows[queryHits(hits)],
  rep[subjectHits(hits)]
)

cov_df <- data.frame(
  win = queryHits(hits),
  bp = width(ov)
)

bp_cov <- tapply(cov_df$bp, cov_df$win, sum)

windows$repeat_bp <- 0
windows$repeat_bp[as.integer(names(bp_cov))] <- bp_cov
windows$repeat_pct <- windows$repeat_bp / width(windows) * 100

# -----------------------------
# Mark windows overlapping SV regions
# -----------------------------

sv_gr <- GRanges(
  seqnames = sv_regions$chr,
  ranges = IRanges(start = sv_regions$start, end = sv_regions$end),
  sv_id = sv_regions$sv_id
)

windows$in_sv <- countOverlaps(windows, sv_gr) > 0

df <- as.data.frame(windows)
df$mid <- (df$start + df$end) / 2 / 1e6

# -----------------------------
# Summary statistics
# -----------------------------

summary_df <- df %>%
  group_by(seqnames, in_sv) %>%
  summarise(
    n_windows = n(),
    median_repeat = median(repeat_pct, na.rm = TRUE),
    mean_repeat = mean(repeat_pct, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  summary_df,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Wilcoxon tests
# -----------------------------

wilcox_df <- data.frame()

for (chr in names(chrom_lengths)) {
  
  df_chr <- subset(df, seqnames == chr)
  test <- wilcox.test(repeat_pct ~ in_sv, data = df_chr)
  
  wilcox_df <- rbind(
    wilcox_df,
    data.frame(
      chromosome = chr,
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  )
}

write.table(
  wilcox_df,
  file = wilcox_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Plot repeat density along chromosomes
# -----------------------------

pdf(plot_output, width = 10, height = 5)

for (chr in names(chrom_lengths)) {
  
  df_chr <- subset(df, seqnames == chr)
  sv_chr <- subset(sv_regions, chr == !!chr)
  
  p <- ggplot(df_chr, aes(x = mid, y = repeat_pct)) +
    geom_line() +
    geom_rect(
      data = sv_chr,
      aes(
        xmin = start / 1e6,
        xmax = end / 1e6,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    theme_minimal() +
    labs(
      title = chr,
      x = "Position (Mb)",
      y = "Repeat density (%)"
    )
  
  print(p)
}

dev.off()
