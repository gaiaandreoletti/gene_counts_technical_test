# gene_counts_technical_test

## Overview

An RNA-seq quality control and comparative analysis framework for evaluating 4 different lab conditions (A, B, C, D) at a clinical laboratory.

**Research question:** Which of the 4 lab conditions performs best? Are any equally good, or should some be dropped?

![image](https://user-images.githubusercontent.com/46945609/227075619-bfceac21-6338-4bdd-b77d-f200d8752f09.png)

---

## Dataset

| File | Description |
|------|-------------|
| `gene_counts_technical_test.csv` | Raw gene count matrix: 60,625 genes × 24 samples (3.8 MB) |
| `Sample_metrics.csv` | Per-sample QC metadata: 13 metrics × 24 samples |

**Experimental design:** 4 lab conditions × 3 biological replicates × 2 sequencing runs = 24 samples total.

All counts are raw counts from a standard RNA-seq pipeline: trim reads → map reads → remove duplicates → convert to counts.

---

## Tech Stack

- **Language:** R
- **Statistical analysis:** DESeq2, DaMiRseq, rstatix
- **Visualization:** ggplot2, pheatmap, factoextra, raincloudplots, ggpubr, ggrepel, gridExtra
- **Data manipulation:** dplyr, tidyverse, reshape2
- **Reporting:** arsenal

---

## Main Entry Points

| File | Description |
|------|-------------|
| `Gaia_mirvie.R` | Primary analysis script (468 lines) — run this for the full pipeline |
| `For_Mirvie.Rmd` | R Markdown version (345 lines) — renders a full HTML/PDF report |

---

## Analysis Workflow

1. Load count matrix and QC metadata
2. Compute summary statistics per condition
3. Visualize 13 QC metrics (read counts, bias metrics, base-type percentages)
4. DESeq2 normalization → PCA and hierarchical clustering
5. Gene filtering by CPM thresholds (0, 10, 50)
6. Heatmap of the 1,000 most variable genes
7. Cook's distance outlier diagnostics
8. DaMiRseq QC pipeline
9. Statistical comparisons with p-values across conditions

---

## Outputs

Running the scripts generates the following PDF reports:

| Output | Description |
|--------|-------------|
| `PCA_QC.pdf` | PCA plots coloured by condition, run, and replicate |
| `heatmap_1000MostVariableGenes.pdf` | Hierarchical clustering of top 1,000 variable genes |
| `genes_cpm_above_*.pdf` | Gene expression filtering at different CPM thresholds |
| `Biases.pdf` | 3′ bias and 5′/3′ bias comparisons across conditions |
| `DAMIRseq_qc_plots.pdf` | DaMiRseq QC pipeline outputs |
| `abcc2_plot*.pdf` | ABCC2 gene-specific analysis |

---

## Author

Gaia Andreoletti — April 2022
