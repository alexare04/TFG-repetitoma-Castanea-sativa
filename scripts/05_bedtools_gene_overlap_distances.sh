#!/bin/bash

# ============================================================
# Script: 05_bedtools_gene_overlap_distances.sh
# Purpose: Calculate physical overlap and distance between
#          repetitive elements / LTR retrotransposons and genes.
# Author: Alejandro Arévalo Sánchez
# ============================================================

set -euo pipefail

# -----------------------------
# Input files
# -----------------------------

GENES_GFF="data/genomic.gff"
LTR_INTACT_GFF="data/castano_genoma_filtrado.fasta.pass.list.gff3"
LTR_TOTAL_GFF="data/castano_genoma_filtrado.fasta.out.gff3"
REPEATMASKER_GFF="data/castano_genoma_filtrado.fasta.mod.out.gff"

# -----------------------------
# Output directory
# -----------------------------

OUTDIR="results/gene_overlap_distances"
mkdir -p "$OUTDIR"

# ============================================================
# 1. Convert gene annotation to BED
# ============================================================

grep -P "\tgene\t" "$GENES_GFF" | \
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' \
> "$OUTDIR/genes.bed"

sort -k1,1 -k2,2n "$OUTDIR/genes.bed" \
> "$OUTDIR/genes.sorted.bed"

GENES="$OUTDIR/genes.sorted.bed"

# ============================================================
# 2. Intact LTR retrotransposons: Gypsy and Copia
# ============================================================

awk '$3=="Gypsy_LTR_retrotransposon" {OFS="\t"; print $1, $4-1, $5}' "$LTR_INTACT_GFF" \
> "$OUTDIR/gypsy_intact.bed"

awk '$3=="Copia_LTR_retrotransposon" {OFS="\t"; print $1, $4-1, $5}' "$LTR_INTACT_GFF" \
> "$OUTDIR/copia_intact.bed"

sort -k1,1 -k2,2n "$OUTDIR/gypsy_intact.bed" > "$OUTDIR/gypsy_intact.sorted.bed"
sort -k1,1 -k2,2n "$OUTDIR/copia_intact.bed" > "$OUTDIR/copia_intact.sorted.bed"

# Physical overlap with genes
bedtools intersect -a "$OUTDIR/gypsy_intact.sorted.bed" -b "$GENES" -wa -u \
> "$OUTDIR/gypsy_intact_overlap_genes.bed"

bedtools intersect -a "$OUTDIR/copia_intact.sorted.bed" -b "$GENES" -wa -u \
> "$OUTDIR/copia_intact_overlap_genes.bed"

# Distance to nearest gene
bedtools closest -a "$OUTDIR/gypsy_intact.sorted.bed" -b "$GENES" -d \
> "$OUTDIR/gypsy_intact_distance.txt"

bedtools closest -a "$OUTDIR/copia_intact.sorted.bed" -b "$GENES" -d \
> "$OUTDIR/copia_intact_distance.txt"

# ============================================================
# 3. Total LTR retrotransposons: Gypsy and Copia
# ============================================================

awk '$0 !~ /^#/ && $0 ~ /Gypsy/ {OFS="\t"; print $1, $4-1, $5}' "$LTR_TOTAL_GFF" \
> "$OUTDIR/gypsy_total.bed"

awk '$0 !~ /^#/ && $0 ~ /Copia/ {OFS="\t"; print $1, $4-1, $5}' "$LTR_TOTAL_GFF" \
> "$OUTDIR/copia_total.bed"

sort -k1,1 -k2,2n "$OUTDIR/gypsy_total.bed" > "$OUTDIR/gypsy_total.sorted.bed"
sort -k1,1 -k2,2n "$OUTDIR/copia_total.bed" > "$OUTDIR/copia_total.sorted.bed"

# Physical overlap with genes
bedtools intersect -a "$OUTDIR/gypsy_total.sorted.bed" -b "$GENES" -wa -u \
> "$OUTDIR/gypsy_total_overlap_genes.bed"

bedtools intersect -a "$OUTDIR/copia_total.sorted.bed" -b "$GENES" -wa -u \
> "$OUTDIR/copia_total_overlap_genes.bed"

# Distance to nearest gene
bedtools closest -a "$OUTDIR/gypsy_total.sorted.bed" -b "$GENES" -d \
> "$OUTDIR/gypsy_total_distance.txt"

bedtools closest -a "$OUTDIR/copia_total.sorted.bed" -b "$GENES" -d \
> "$OUTDIR/copia_total_distance.txt"

# ============================================================
# 4. All repetitive elements from RepeatMasker
# ============================================================

awk '$0 !~ /^#/ {OFS="\t"; print $1, $4-1, $5}' "$REPEATMASKER_GFF" \
> "$OUTDIR/repeats_total.bed"

sort -k1,1 -k2,2n "$OUTDIR/repeats_total.bed" \
> "$OUTDIR/repeats_total.sorted.bed"

# Genes with at least one physically overlapping repetitive element
bedtools intersect -a "$GENES" -b "$OUTDIR/repeats_total.sorted.bed" -wa -u \
> "$OUTDIR/genes_with_repeats.bed"

# Fraction of each gene covered by repetitive elements
bedtools coverage -a "$GENES" -b "$OUTDIR/repeats_total.sorted.bed" \
> "$OUTDIR/gene_repeat_coverage.txt"

# Distance from each repetitive element to nearest gene
bedtools closest -a "$OUTDIR/repeats_total.sorted.bed" -b "$GENES" -d \
> "$OUTDIR/repeats_distance.txt"

# ============================================================
# 5. Summary counts
# ============================================================

SUMMARY="$OUTDIR/overlap_summary.tsv"

echo -e "category\ttotal_elements\toverlapping_genes_elements" > "$SUMMARY"

echo -e "Gypsy_intact\t$(wc -l < "$OUTDIR/gypsy_intact.sorted.bed")\t$(wc -l < "$OUTDIR/gypsy_intact_overlap_genes.bed")" >> "$SUMMARY"
echo -e "Copia_intact\t$(wc -l < "$OUTDIR/copia_intact.sorted.bed")\t$(wc -l < "$OUTDIR/copia_intact_overlap_genes.bed")" >> "$SUMMARY"
echo -e "Gypsy_total\t$(wc -l < "$OUTDIR/gypsy_total.sorted.bed")\t$(wc -l < "$OUTDIR/gypsy_total_overlap_genes.bed")" >> "$SUMMARY"
echo -e "Copia_total\t$(wc -l < "$OUTDIR/copia_total.sorted.bed")\t$(wc -l < "$OUTDIR/copia_total_overlap_genes.bed")" >> "$SUMMARY"

echo "Analysis completed. Results saved in $OUTDIR"
