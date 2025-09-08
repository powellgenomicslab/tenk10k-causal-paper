# 2-mr: Mendelian Randomization Analysis

This directory contains scripts for Mendelian Randomization (MR) analyses and visualization of results comparing TenK10K single-cell eQTL data with bulk eQTL data.

## Overview

This section performs comprehensive Mendelian Randomization analyses to identify causal relationships between gene expression and complex traits. The analysis compares results from TenK10K single-cell eQTL data with eQTLGen bulk tissue results to understand cell-type-specific causal effects.

## Contents

### Scripts
- `mr_results_main.R` - Main figures and analyses for MR results:
  - Creates Figure 2 showing MR overview across cell types and traits
  - Generates tally plots comparing MR-only vs MR & GWAS associations
  - Compares TenK10K single-cell results with eQTLGen bulk results
  - Creates comprehensive visualization panels for manuscript

- `mr_results_supp.R` - Supplementary MR analyses:
  - Additional statistical comparisons
  - Extended visualizations for supplementary materials
  - Detailed breakdowns by cell type and trait categories

## Dependencies

### R packages
```r
library(patchwork)      # Combining multiple plots
library(ragg)           # High-quality graphics device
library(scales)         # Scale functions for visualization
library(paletteer)      # Color palettes
library(ggnewscale)     # Multiple color scales in ggplot2
library(geomtextpath)   # Text along paths in plots
# Plus all dependencies from 0-preprocess/preprocess_results.R
```

## Input Files

Both scripts source `scripts/0-preprocess/preprocess_results.R`, which provides:
- Filtered mSMR results from TenK10K and eQTLGen
- Cell type and trait mappings
- Gene annotations and categories

## Key Analyses

### Main Results (`mr_results_main.R`)
1. **Panel A1** - Tally by number of cell types and GWAS associations
2. **Panel A2** - Tally by number of cell types and eQTLGen MR results
3. **Comparison plots** - TenK10K vs eQTLGen MR effects
4. **Cell-type specificity** - Analysis of cell-type-specific causal effects

### Supplementary Results (`mr_results_supp.R`)
1. Extended statistical comparisons
2. Additional visualization panels
3. Detailed breakdowns by trait categories
4. Quality control and sensitivity analyses

## Output

Generates:
- **Figure 2** - Main MR overview figure for manuscript
- **Supplementary figures** - Extended MR analyses
- **Statistical summaries** - MR effect comparisons
- **Cell-type specificity metrics** - Quantification of cell-type-specific effects

## Usage

Run the MR analysis scripts from the repository root:

```bash
# Main MR results and Figure 2
Rscript scripts/2-mr/mr_results_main.R

# Supplementary MR analyses
Rscript scripts/2-mr/mr_results_supp.R
```

**Prerequisites:** 
1. Complete the Snakemake workflow to generate mSMR results
2. Run the preprocessing script (`0-preprocess/preprocess_results.R`)

## Key Features

- **Multi-dataset comparison** - TenK10K single-cell vs eQTLGen bulk
- **Cell-type resolution** - Analysis at individual cell type level
- **Effect size comparison** - Statistical comparison of MR effect sizes
- **Visualization** - Comprehensive plots for manuscript figures

## Color Schemes

The scripts use consistent color schemes:
- `mr_only` - Light blue (#A4ABB0FF) for MR-only associations
- `mr_other` - Dark blue (#4C6C94FF) for MR & GWAS associations