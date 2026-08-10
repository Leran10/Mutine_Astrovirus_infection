# Murine Astrovirus Infection - scRNA-seq Analysis

This repository contains the analysis code for single-cell RNA sequencing (scRNA-seq) data from murine astrovirus infection studies. The scripts support two separate publications as described below.

---

## Paper 1: Astrovirus scRNA-seq Gene Expression Analysis

**Script:** `Astrovirus_scRNA_gene_expression_final_plots_github.Rmd`

This R Markdown notebook contains the full scRNA-seq analysis pipeline for characterizing the cellular response to murine astrovirus infection in intestinal tissue. The analysis covers data preprocessing, cell clustering, cell-type annotation (using SingleR with ImmGen and SCSA), separation of Epcam+ (epithelial) and Ptprc+/CD45+ (immune) compartments, gene expression visualization, and differential expression analysis across experimental conditions:

- Wild-type naive (WT naive)
- Wild-type fecal filtrate transfer (WT + FFT)
- Rag2⁻/⁻Il2rg⁻/⁻
- Rag2⁻/⁻Il2rg⁻/⁻Ifnlr1⁻/⁻

### Cell Types Analyzed

- **Epithelial cells**: Enterocytes (mature/immature, proximal/distal), Goblet cells, Paneth cells, Tuft cells, Enteroendocrine cells
- **Progenitor cells**: Transit-amplifying/Stem cells, Enterocyte progenitors
- **Immune cells**: T cells, B cells, Macrophages
- **Other**: Stromal cells, Mesenchymal cells, Erythroid cells

---

## Paper 2: Host IFN-λ signaling is required for CD8⁺ T cell-mediated clearance of chronic murine astrovirus infection

**Scripts:** `T cell update CD45 only/`

This directory contains the T cell subtype analysis restricted to the CD45+ immune compartment. T cells were subset from CD45+ clusters only (excluding T cells in Epcam+ cluster 20) and re-analyzed independently using SCTransform normalization, PCA, tSNE, and UMAP dimensionality reduction. Unsupervised re-clustering at resolution 0.5 was followed by T cell subtype assignment based on Cd3e, Cd4, Cd8a, and Cd8b1 marker expression thresholds, yielding CD4+, CD8αα, CD8αβ, mixed CD4/CD8, and double-negative (DN) T cell populations. The analysis compares WT naive and WT + FFT conditions (6,429 T cells total).

| Script | Description |
|--------|-------------|
| `run_full_analysis_CD45_Tcells.R` | Main pipeline: normalization, clustering, subtype labeling, DE analysis, volcano plots, and pathway analysis (ORA + GSEA) |
| `run_remaining_plots.R` | Supplementary visualizations: violin plots, dot plots, feature plots by subtype |
| `run_remaining_plots_2.R` | Additional plots: resolution comparisons, per-cluster visualizations |
| `run_bubble_plot.R` | Bubble plots of marker expression across subtypes |
| `run_fix_plots.R` | Plot refinements and formatting corrections |

---

## Requirements

### R Packages
- Seurat (v4+)
- SeuratDisk
- ggplot2
- ggpubr
- ggrepel
- tidyverse / dplyr / data.table
- fgsea
- msigdbr
- clusterProfiler
- org.Mm.eg.db
- viridis
- svglite

## Data

The analysis requires preprocessed Seurat objects (`.h5seurat` files). These files should be placed in the appropriate directory as specified in the scripts.

## Author

Leran Wang

## License

Please contact the author for licensing information.
