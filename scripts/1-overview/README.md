# 1-overview: Study Overview and Design

This directory contains scripts for generating overview statistics and visualizations presented in the manuscript

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

## Usage

Run the overview script from the repository root:

```bash
Rscript scripts/1-overview/study_design.R
```
