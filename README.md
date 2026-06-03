# Computational analysis of repetitive DNA in *Castanea sativa*

This repository contains the scripts, summary tables and figures generated during the Bachelor's Thesis:

**Computational analysis of repetitive DNA in the European chestnut (*Castanea sativa* Mill.)**

The main aim of this work was to characterize the repetitive DNA fraction of the *Castanea sativa* reference genome, with special focus on LTR retrotransposons from the *Gypsy* and *Copia* superfamilies.

## Repository structure

```text
.
├── scripts/
│   ├── 00_run_repeatmodeler_repeatmasker.sh
│   ├── 01_bedtools_gene_overlap_distances.sh
│   ├── 02_plot_repeat_distance_to_genes.R
│   ├── 03_structural_variants_repeat_density.R
│   ├── 04_run_edta_ltr_retriever.sh
│   ├── 05_ltr_distribution_windows.R
│   ├── 06_gene_ltr_correlation.R
│   ├── 07_go_enrichment_ltr_genes.R
│   ├── 08_ltr_insertion_age.py
│   └── 09_plot_ltr_insertion_age.R
│
├── results/
│   ├── figures/
│   └── tables/
│
├── environment.yml
└── README.md
```

## Scripts

The scripts included in this repository reproduce the main analyses performed in the thesis.

### `00_run_repeatmodeler_repeatmasker.sh`

Runs RepeatModeler and RepeatMasker to generate a repeat library and obtain a genome-wide annotation of repetitive DNA.

### `01_bedtools_gene_overlap_distances.sh`

Calculates the physical overlap between genes and repetitive elements using `bedtools intersect`, as well as the distance from repetitive elements to the nearest annotated gene using `bedtools closest`. It also groups RepeatMasker annotations by repeat class and calculates the proportion of genes affected by each repeat group.

### `02_plot_repeat_distance_to_genes.R`

Plots the distance from different classes of repetitive elements to the nearest gene, considering only non-overlapping repetitive elements.

### `03_structural_variants_repeat_density.R`

Analyzes repetitive element density inside and outside large structural variant regions detected between haplotypes.

### `04_run_edta_ltr_retriever.sh`

Runs EDTA for transposable element annotation and identification of intact LTR retrotransposons using LTR_retriever.

### `05_ltr_distribution_windows.R`

Calculates the distribution of *Gypsy* and *Copia* LTR retrotransposons and gene annotations in 1 Mb windows along the chromosomes.

### `06_gene_ltr_correlation.R`

Calculates Spearman correlations between gene density and the abundance of *Gypsy* and *Copia* LTR retrotransposons.

### `07_go_enrichment_ltr_genes.R`

Performs Gene Ontology enrichment analysis for genes associated with LTR retrotransposons, separately for genes associated with *Gypsy* and *Copia* elements.

### `08_ltr_insertion_age.py`

Estimates the insertion age of intact LTR retrotransposons from the identity value between their two LTR regions.

### `09_plot_ltr_insertion_age.R`

Plots insertion age distributions of intact LTR retrotransposons located inside and outside genes and calculates summary statistics by superfamily and location.

## Input files

Large input files are not included in this repository due to file size limitations. To run the scripts, the required genome assembly, gene annotation and program output files should be placed in a `data/` directory.

Expected main input files:

```text
data/
├── castano_genoma_filtrado.fasta
├── genomic.gff
├── castano_genoma_filtrado.fasta.out.gff3
├── castano_genoma_filtrado.fasta.pass.list.gff3
├── castano_genoma_filtrado.fasta.LTRlib.fa
├── castano_genoma_filtrado.fasta.mod.out
├── castano_genoma_filtrado.fasta.mod.out.gff
├── gene2go.tsv
├── all_genes_mappable.txt
├── genes_with_Gypsy_mappable.txt
└── genes_with_Copia_mappable.txt
```

## Results

The `results/` directory contains the summary tables and final figures generated during the analysis.

### Main figures

```text
results/figures/
├── GO_enrichment_Gypsy_Copia.pdf
├── Gypsy_Copia_gene_density_correlations.pdf
├── LTR_Gypsy_Copia_Genes_by_chromosome.pdf
├── LTR_insertion_age_inside_outside_genes.pdf
├── SV_repeat_density_by_chromosome.pdf
└── repeat_distance_to_genes_by_class.pdf
```

### Main tables

```text
results/tables/
├── GO_enrichment_Copia.tsv
├── GO_enrichment_Gypsy.tsv
├── SV_repeat_density_summary.tsv
├── SV_repeat_density_wilcoxon.tsv
├── chromosome_spearman_gene_ltr.tsv
├── global_spearman_gene_ltr.tsv
├── ltr_insertion_age_summary.tsv
├── repeat_distance_to_genes_summary.tsv
└── window_counts_1Mb.tsv
```

## Software environment

The main tools and packages used in this work are listed in `environment.yml`.

To create the conda environment:

```bash
conda env create -f environment.yml
conda activate castanea_repetitive_dna
```

Some tools, such as EDTA, RepeatMasker and RepeatModeler, may require additional configuration depending on the operating system and local installation.

## Reproducibility note

This repository contains the main scripts, summary tables and final figures used in the thesis. Complete genome files and large intermediate files are not included, but they can be regenerated by running the scripts in the indicated order with the corresponding input files.

## Author

Alejandro Arévalo Sánchez
Bachelor's Thesis
Degree in Biotechnology
