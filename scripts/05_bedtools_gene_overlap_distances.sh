#!/bin/bash

# ============================================================
# Script: 05_bedtools_gene_overlap_distances.sh
# Purpose: Calculate physical overlap and distance between
#          repetitive elements / LTR retrotransposons and genes.
#          It also groups RepeatMasker annotations by repeat class
#          and calculates gene overlap and nearest-gene distance
#          for each repeat group.
# Author: Alejandro Arévalo Sánchez
# ============================================================

set -euo pipefail

# -----------------------------
# Input files
# -----------------------------

GENES_GFF="data/genomic.gff"

LTR_INTACT_GFF="data/castano_genoma_filtrado.fasta.pass.list.gff3"
LTR_TOTAL_GFF="data/castano_genoma_filtrado.fasta.out.gff3"

# RepeatMasker GFF-like annotation used for global repeat overlap
REPEATMASKER_GFF="data/castano_genoma_filtrado.fasta.mod.out.gff"

# RepeatMasker .out file used to recover repeat classes
REPEATMASKER_OUT="data/castano_genoma_filtrado.fasta.mod.out"

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

sort -k1,1 -k2,2n "$OUTDIR/gypsy_intact.bed" \
> "$OUTDIR/gypsy_intact.sorted.bed"

sort -k1,1 -k2,2n "$OUTDIR/copia_intact.bed" \
> "$OUTDIR/copia_intact.sorted.bed"

# Physical overlap with genes
bedtools intersect \
  -a "$OUTDIR/gypsy_intact.sorted.bed" \
  -b "$GENES" \
  -wa -u \
> "$OUTDIR/gypsy_intact_overlap_genes.bed"

bedtools intersect \
  -a "$OUTDIR/copia_intact.sorted.bed" \
  -b "$GENES" \
  -wa -u \
> "$OUTDIR/copia_intact_overlap_genes.bed"

# Distance to nearest gene
bedtools closest \
  -a "$OUTDIR/gypsy_intact.sorted.bed" \
  -b "$GENES" \
  -d \
> "$OUTDIR/gypsy_intact_distance.txt"

bedtools closest \
  -a "$OUTDIR/copia_intact.sorted.bed" \
  -b "$GENES" \
  -d \
> "$OUTDIR/copia_intact_distance.txt"

# ============================================================
# 3. Total LTR retrotransposons: Gypsy and Copia
# ============================================================

awk '$0 !~ /^#/ && $0 ~ /Gypsy/ {OFS="\t"; print $1, $4-1, $5}' "$LTR_TOTAL_GFF" \
> "$OUTDIR/gypsy_total.bed"

awk '$0 !~ /^#/ && $0 ~ /Copia/ {OFS="\t"; print $1, $4-1, $5}' "$LTR_TOTAL_GFF" \
> "$OUTDIR/copia_total.bed"

sort -k1,1 -k2,2n "$OUTDIR/gypsy_total.bed" \
> "$OUTDIR/gypsy_total.sorted.bed"

sort -k1,1 -k2,2n "$OUTDIR/copia_total.bed" \
> "$OUTDIR/copia_total.sorted.bed"

# Physical overlap with genes
bedtools intersect \
  -a "$OUTDIR/gypsy_total.sorted.bed" \
  -b "$GENES" \
  -wa -u \
> "$OUTDIR/gypsy_total_overlap_genes.bed"

bedtools intersect \
  -a "$OUTDIR/copia_total.sorted.bed" \
  -b "$GENES" \
  -wa -u \
> "$OUTDIR/copia_total_overlap_genes.bed"

# Distance to nearest gene
bedtools closest \
  -a "$OUTDIR/gypsy_total.sorted.bed" \
  -b "$GENES" \
  -d \
> "$OUTDIR/gypsy_total_distance.txt"

bedtools closest \
  -a "$OUTDIR/copia_total.sorted.bed" \
  -b "$GENES" \
  -d \
> "$OUTDIR/copia_total_distance.txt"

# ============================================================
# 4. All repetitive elements from RepeatMasker
# ============================================================

awk '$0 !~ /^#/ {OFS="\t"; print $1, $4-1, $5}' "$REPEATMASKER_GFF" \
> "$OUTDIR/repeats_total.bed"

sort -k1,1 -k2,2n "$OUTDIR/repeats_total.bed" \
> "$OUTDIR/repeats_total.sorted.bed"

# Genes with at least one physically overlapping repetitive element
bedtools intersect \
  -a "$GENES" \
  -b "$OUTDIR/repeats_total.sorted.bed" \
  -wa -u \
> "$OUTDIR/genes_with_repeats.bed"

# Fraction of each gene covered by repetitive elements
bedtools coverage \
  -a "$GENES" \
  -b "$OUTDIR/repeats_total.sorted.bed" \
> "$OUTDIR/gene_repeat_coverage.txt"

# Distance from each repetitive element to nearest gene
bedtools closest \
  -a "$OUTDIR/repeats_total.sorted.bed" \
  -b "$GENES" \
  -d \
> "$OUTDIR/repeats_distance.txt"

# ============================================================
# 5. Summary counts for LTR overlap
# ============================================================

SUMMARY="$OUTDIR/overlap_summary.tsv"

echo -e "category\ttotal_elements\toverlapping_gene_elements\tpct_overlapping_gene_elements" > "$SUMMARY"

for category in gypsy_intact copia_intact gypsy_total copia_total; do

  total=$(wc -l < "$OUTDIR/${category}.sorted.bed")
  overlap=$(wc -l < "$OUTDIR/${category}_overlap_genes.bed")

  pct=$(awk -v overlap="$overlap" -v total="$total" \
    'BEGIN{
      if(total > 0) printf "%.4f", (overlap/total)*100;
      else printf "NA"
    }')

  echo -e "${category}\t${total}\t${overlap}\t${pct}" >> "$SUMMARY"

done

# ============================================================
# 6. Global repeat overlap summary
# ============================================================

GENE_REPEAT_SUMMARY="$OUTDIR/gene_repeat_coverage_summary.tsv"

total_genes=$(wc -l < "$GENES")
genes_with_repeats=$(wc -l < "$OUTDIR/genes_with_repeats.bed")

pct_genes_with_repeats=$(awk -v n="$genes_with_repeats" -v N="$total_genes" \
  'BEGIN{
    if(N > 0) printf "%.4f", (n/N)*100;
    else printf "NA"
  }')

# In bedtools coverage output, column 7 is the fraction of each A interval covered by B.
mean_gene_fraction=$(awk '{sum+=$7; n++} END{if(n>0) printf "%.6f", sum/n; else print "NA"}' \
  "$OUTDIR/gene_repeat_coverage.txt")

mean_gene_pct=$(awk -v f="$mean_gene_fraction" \
  'BEGIN{
    if(f!="NA") printf "%.4f", f*100;
    else printf "NA"
  }')

echo -e "metric\tvalue" > "$GENE_REPEAT_SUMMARY"
echo -e "total_genes\t${total_genes}" >> "$GENE_REPEAT_SUMMARY"
echo -e "genes_with_repeats\t${genes_with_repeats}" >> "$GENE_REPEAT_SUMMARY"
echo -e "pct_genes_with_repeats\t${pct_genes_with_repeats}" >> "$GENE_REPEAT_SUMMARY"
echo -e "mean_gene_fraction_covered_by_repeats\t${mean_gene_fraction}" >> "$GENE_REPEAT_SUMMARY"
echo -e "mean_gene_pct_covered_by_repeats\t${mean_gene_pct}" >> "$GENE_REPEAT_SUMMARY"

# ============================================================
# 7. Repeat classes: genes affected and distance to nearest gene
# ============================================================

GROUPDIR="$OUTDIR/RM_geneImpact_grouped"
mkdir -p "$GROUPDIR"

# -----------------------------
# Convert RepeatMasker .out file to BED with repeat class
# -----------------------------

awk 'BEGIN{OFS="\t"}
     NR<=3 {next}
     $1=="SW" || $1=="score" || $1=="" {next}
     {
       chr=$5
       start=$6-1
       end=$7
       strand=$9
       rep=$10
       class=$11
       print chr, start, end, class, rep, strand
     }' "$REPEATMASKER_OUT" \
> "$GROUPDIR/repeats.byclass.bed"

sort -k1,1 -k2,2n "$GROUPDIR/repeats.byclass.bed" \
> "$GROUPDIR/repeats.byclass.sorted.bed"

# -----------------------------
# Group repeat classes into broader categories
# -----------------------------

awk 'BEGIN{OFS="\t"}
{
  cls=$4
  grp="OTHER"

  if(cls=="LTR/Gypsy") grp="Gypsy"
  else if(cls=="LTR/Copia") grp="Copia"
  else if(cls ~ /^LINE\//) grp="LINE"
  else if(cls ~ /^DNA/) grp="DNA_transposons"
  else if(cls=="RC/Helitron") grp="Rolling_circles"
  else if(cls=="Simple_repeat") grp="Simple_repeats"
  else if(cls=="Low_complexity") grp="Low_complexity"
  else if(cls=="tRNA" || cls=="rRNA" || cls=="snRNA") grp="Small_RNA"
  else if(cls=="unknown" || cls=="LTR/unknown") grp="Unknown"
  else next

  print $1,$2,$3,grp,$5,$6
}' "$GROUPDIR/repeats.byclass.sorted.bed" \
> "$GROUPDIR/repeats.grouped.bed"

sort -k1,1 -k2,2n "$GROUPDIR/repeats.grouped.bed" \
> "$GROUPDIR/repeats.grouped.sorted.bed"

cut -f4 "$GROUPDIR/repeats.grouped.sorted.bed" | sort -u \
> "$GROUPDIR/groups.txt"

# -----------------------------
# Percentage of genes affected by each repeat group
# -----------------------------

NGENES=$(wc -l < "$GENES")

echo -e "group\tgenes_affected\tpct_genes_affected" \
> "$GROUPDIR/genes_affected_by_group.tsv"

while read -r grp; do

  awk -v G="$grp" '$4==G' "$GROUPDIR/repeats.grouped.sorted.bed" \
  > "$GROUPDIR/${grp}.bed"

  n=$(bedtools intersect \
        -a "$GENES" \
        -b "$GROUPDIR/${grp}.bed" \
        -wa -u | wc -l)

  pct=$(awk -v n="$n" -v N="$NGENES" \
    'BEGIN{
      if(N > 0) printf "%.4f", (n/N)*100;
      else printf "NA"
    }')

  echo -e "${grp}\t${n}\t${pct}" \
  >> "$GROUPDIR/genes_affected_by_group.tsv"

done < "$GROUPDIR/groups.txt"

# -----------------------------
# Distance to nearest gene by repeat group
# -----------------------------

echo -e "group\tdistance_bp" \
> "$GROUPDIR/distances_by_group.long.tsv"

while read -r grp; do

  bedtools closest \
    -a "$GROUPDIR/${grp}.bed" \
    -b "$GENES" \
    -d \
  > "$GROUPDIR/closest_${grp}.txt"

  awk -v G="$grp" 'BEGIN{OFS="\t"} {print G, $NF}' \
  "$GROUPDIR/closest_${grp}.txt" \
  >> "$GROUPDIR/distances_by_group.long.tsv"

done < "$GROUPDIR/groups.txt"

# -----------------------------
# Summary of distances by repeat group
# -----------------------------

echo -e "group\tn_repeats\tmin\tq1\tmedian\tmean\tq3\tmax" \
> "$GROUPDIR/distances_by_group.summary.tsv"

while read -r grp; do

  gawk -v G="$grp" '
  {d[NR]=$NF; sum+=$NF}
  END{
    if(NR==0){exit}
    n=NR
    asort(d)
    min=d[1]
    max=d[n]
    mean=sum/n
    q1=d[int((n*0.25)<1?1:(n*0.25))]
    med=d[int((n*0.50)<1?1:(n*0.50))]
    q3=d[int((n*0.75)<1?1:(n*0.75))]
    printf "%s\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\n", G, n, min, q1, med, mean, q3, max
  }' "$GROUPDIR/closest_${grp}.txt" \
  >> "$GROUPDIR/distances_by_group.summary.tsv"

done < "$GROUPDIR/groups.txt"

# ============================================================
# 8. Final message
# ============================================================

echo "Analysis completed."
echo "Main output directory: $OUTDIR"
echo "Repeat class output directory: $GROUPDIR"
