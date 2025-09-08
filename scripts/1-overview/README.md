# 1-overview: Study Overview and Design

This directory contains scripts for generating overview statistics and visualizations for the TenK10K causal inference study.

## Overview

This section provides summary statistics and creates figures that give an overview of the study design, data characteristics, and analytical approach. It generates key numbers and visualizations used in the manuscript to describe the scope and scale of the analyses.

## Contents

### Scripts
- `study_design.R` - Main script that:
  - Sources preprocessed results from `0-preprocess/`
  - Calculates gene counts by cell type and gene category
  - Generates summary statistics for the manuscript
  - Creates visualizations of study design and data characteristics
  - Produces overview figures showing gene distributions across cell types

## Dependencies

### R packages
```r
library(patchwork)      # Combining multiple plots
library(paletteer)      # Color palettes
library(geomtextpath)   # Text along paths in ggplot2
# Plus all dependencies from 0-preprocess/preprocess_results.R:
# data.table, tidyverse, arrow, fs, qvalue
```

## Input Files

The script sources `scripts/0-preprocess/preprocess_results.R`, which loads:
- Preprocessed MAGMA and mSMR results
- Cell type and trait mappings
- Gene annotations

## Output

Generates:
- Summary statistics for manuscript text
- Overview figures for the study design
- Gene count summaries by cell type and gene category (protein-coding vs non-coding)
- Data characteristics tables

## Usage

Run the overview script from the repository root:

```bash
Rscript scripts/1-overview/study_design.R
```

**Prerequisites:** 
1. Complete the Snakemake workflow to generate input data
2. Run the preprocessing script (`0-preprocess/preprocess_results.R`)

## Key Outputs

The script categorizes genes into:
- **Protein-coding genes** - Standard protein-coding loci
- **Non-coding genes** - Non-protein-coding loci (e.g., lncRNA, miRNA)

And provides counts across:
- Individual cell types from TenK10K Phase 1
- Major cell type categories
- Overall gene universe used in the study

## Figures Generated

- Study design overview plots
- Gene distribution visualizations by cell type
- Summary statistics plots for data characteristics