# Murine Astrovirus Infection - scRNA-seq Analysis

This repository contains the analysis code for single-cell RNA sequencing (scRNA-seq) data from murine astrovirus infection studies.

## Contents

- `Astrovirus_scRNA_gene_expression_final_plots_github.Rmd` - R Markdown notebook for scRNA-seq gene expression analysis and visualization

## Overview

This analysis examines the cellular response to astrovirus infection in murine intestinal tissue using single-cell RNA sequencing. The study compares multiple experimental conditions:

- Wild-type naive (Wtnaive)
- Wild-type fecal filtrate transfer (WTFFT)
- Rag2Il2rg knockout
- Rag2Il2rgIfnlr1 knockout

## Cell Types Analyzed

The analysis includes multiple intestinal epithelial and immune cell types:
- **Epithelial cells**: Enterocytes (mature/immature, proximal/distal), Goblet cells, Paneth cells, Tuft cells, Enteroendocrine cells
- **Progenitor cells**: Transit-amplifying/Stem cells, Enterocyte progenitors
- **Immune cells**: T cells, B cells, Macrophages
- **Other**: Stromal cells, Mesenchymal cells, Erythroid cells

## Requirements

### R Packages
- Seurat (v4+)
- SeuratDisk
- ggplot2
- ggpubr
- plotly
- ggrepel
- gridExtra
- viridis
- svglite
- fgsea
- msigdbr
- tidyverse
- dplyr
- data.table
- shiny

## Usage

Open the R Markdown file in RStudio and run the chunks sequentially. The notebook is set up as an interactive Shiny document.

```r
# Open in RStudio
file.edit("Astrovirus_scRNA_gene_expression_final_plots_github.Rmd")
```

## Data

The analysis requires preprocessed Seurat objects (`.h5seurat` files). These files should be placed in the appropriate directory as specified in the script.

## Author

Leran Wang

## Date

August 11, 2022

## License

Please contact the author for licensing information.
