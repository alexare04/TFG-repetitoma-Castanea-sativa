#!/usr/bin/env python3

# ============================================================
# Script: 07_ltr_insertion_age.py
# Purpose: Estimate insertion age of intact LTR retrotransposons
#          located inside and outside annotated genes.
# Author: Alejandro Arévalo Sánchez
# ============================================================

import re
from pathlib import Path

# -----------------------------
# Parameters
# -----------------------------

# Substitution rate from Juglans regia:
# substitutions/site/year
MU = 2.29e-9

# -----------------------------
# Input and output files
# -----------------------------

input_files = {
    "in_genes": Path("results/gene_overlap_distances/intact_in_genes.gff3"),
    "outside_genes": Path("results/gene_overlap_distances/intact_outside_genes.gff3"),
}

output_dir = Path("results/ltr_insertion_age")
output_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Functions
# -----------------------------

def identity_to_age_mya(identity: float, mu: float = MU) -> float:
    """
    Convert LTR identity into estimated insertion age.

    Formula:
        T = (1 - identity) / (2 * mu)

    Output is converted from years to million years ago (MYA).
    """
    return (1 - identity) / (2 * mu) / 1e6


def extract_ltr_identity(line: str):
    """
    Extract ltr_identity value from a GFF3 line.
    """
    match = re.search(r"ltr_identity=([0-9.]+)", line)
    if match:
        return float(match.group(1))
    return None


def classify_ltr(line: str):
    """
    Classify the LTR retrotransposon as Gypsy or Copia.
    """
    if "Gypsy_LTR_retrotransposon" in line:
        return "Gypsy"
    if "Copia_LTR_retrotransposon" in line:
        return "Copia"
    return None


def process_gff(gff_file: Path, location_label: str):
    """
    Process one GFF3 file and write age estimates for Gypsy and Copia.
    """
    output_files = {
        "Gypsy": output_dir / f"Gypsy_{location_label}_age_MYA.txt",
        "Copia": output_dir / f"Copia_{location_label}_age_MYA.txt",
    }

    counts = {"Gypsy": 0, "Copia": 0}

    with gff_file.open() as gff, \
         output_files["Gypsy"].open("w") as gypsy_out, \
         output_files["Copia"].open("w") as copia_out:

        writers = {
            "Gypsy": gypsy_out,
            "Copia": copia_out,
        }

        for line in gff:
            if line.startswith("#"):
                continue

            family = classify_ltr(line)
            if family is None:
                continue

            identity = extract_ltr_identity(line)
            if identity is None:
                continue

            age = identity_to_age_mya(identity)
            writers[family].write(f"{age:.4f}\n")
            counts[family] += 1

    return counts


# -----------------------------
# Main analysis
# -----------------------------

for location_label, gff_file in input_files.items():

    if not gff_file.exists():
        raise FileNotFoundError(f"Input file not found: {gff_file}")

    counts = process_gff(gff_file, location_label)

    print(f"Processed file: {gff_file}")
    print(f"  Gypsy elements: {counts['Gypsy']}")
    print(f"  Copia elements: {counts['Copia']}")

print(f"Results saved in: {output_dir}")
