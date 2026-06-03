#!/bin/bash

# ============================================================
# Script: 00_run_edta_ltr_retriever.sh
# Purpose: Run EDTA and LTR_retriever to identify and classify
#          transposable elements and intact LTR retrotransposons
#          in the Castanea sativa reference genome.
# Author: Alejandro Arévalo Sánchez
# ============================================================

set -euo pipefail

# -----------------------------
# Input file
# -----------------------------

GENOME="data/castano_genoma_filtrado.fasta"

# -----------------------------
# Parameters
# -----------------------------

THREADS=16
SPECIES="others"

# -----------------------------
# Output directory
# -----------------------------

OUTDIR="results/EDTA_LTR_retriever"
mkdir -p "$OUTDIR"

# EDTA writes most output files in the same directory as the genome.
# Therefore, copy the genome to the output directory before running EDTA.

cp "$GENOME" "$OUTDIR/"

GENOME_BASENAME=$(basename "$GENOME")
GENOME_IN_OUTDIR="$OUTDIR/$GENOME_BASENAME"

# -----------------------------
# Run EDTA
# -----------------------------

EDTA.pl \
  --genome "$GENOME_IN_OUTDIR" \
  --species "$SPECIES" \
  --threads "$THREADS" \
  --step all

echo "EDTA analysis completed."
echo "Output directory: $OUTDIR"
echo "Check the generated EDTA and LTR_retriever files for TE annotations and intact LTR retrotransposons."
