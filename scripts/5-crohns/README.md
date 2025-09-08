# 5-crohns: Crohn's Disease Case Study

This directory contains scripts for the Crohn's disease case study, which provides a detailed example of applying the causal inference framework to a specific disease with matched single-cell data.

## Overview

This section presents a comprehensive case study of Crohn's disease, integrating genetic associations, single-cell differential expression analysis, and causal inference. The analysis demonstrates how the TenK10K causal inference framework can be applied to understand disease mechanisms at cellular resolution.

## Contents

### Subdirectories

#### `deg/` - Differential Expression Gene Analysis
- `1-prepare_crohns_deg_pre-processed.R` - Prepare Crohn's DEG data for analysis
- `2-find_deg_mr_overlap_genes_both_colon_and_ti.R` - Find overlapping genes between DEG and MR results (colon and terminal ileum)
- `3-prepare_data_for_MR_comparison.py` - Prepare data for MR comparison analysis
- `4-plot_UMAP_cell_annotation.py` - Create UMAP plots with cell type annotations
- `5-plot_UMAP_genes_of_interest_counts.py` - Plot expression of genes of interest on UMAP

#### `figures/` - Figure Generation
- `Crohns_example_locus_zoom.R` - Create locus zoom plots for Crohn's disease examples
- `Figure5-combined_Crohns_figure.R` - Generate combined Figure 5 for manuscript
- `annotated_heatmap_only_canonicalanddeg.R` - Create annotated heatmaps for canonical and DEG results
- `barplot_comparison_with_eqtlgen_and_magma.R` - Comparison bar plots with eQTLGen and MAGMA
- `barplot_gene_numbers_by_celltype.R` - Bar plots showing gene numbers by cell type
- `crohns_deg_mr_comparison.R` - Comprehensive comparison of DEG and MR results

#### `prepare_data/` - Data Preparation
Scripts for preparing Crohn's disease-specific datasets

#### `supplementary/` - Supplementary Analyses
Additional analyses and supplementary figures for Crohn's disease case study

## Key Analyses

### 1. Differential Expression Analysis
- Identification of genes differentially expressed in Crohn's disease
- Cell-type-specific differential expression patterns
- Integration with genetic association signals

### 2. MR-DEG Integration
- Overlap analysis between MR causal genes and DEGs
- Validation of causal effects using expression changes
- Cell-type-specific validation in relevant tissues (colon, terminal ileum)

### 3. Visualization
- UMAP plots showing cell type annotations and gene expression
- Locus zoom plots for key genetic associations
- Heatmaps integrating multiple data types
- Comparison plots with external datasets

### 4. Validation Analysis
- Comparison with eQTLGen bulk tissue results
- Integration with MAGMA gene-level associations
- Cross-validation across different analytical approaches

## Dependencies

### R packages
```r
# Standard packages from preprocessing
# Visualization packages for plots and heatmaps
# Statistical packages for DEG analysis
```

### Python packages
```python
# Single-cell analysis (scanpy, pandas, numpy)
# Visualization (matplotlib, seaborn)
# Data manipulation and analysis
```

## Input Files

The analysis uses:
- Crohn's disease GWAS summary statistics
- Single-cell RNA-seq data from Crohn's disease samples
- TenK10K MR results for Crohn's disease
- Cell type annotations and metadata
- Reference datasets (eQTLGen, MAGMA)

## Output

Generates:
- **Figure 5** - Main Crohn's disease case study figure
- **DEG-MR overlap results** - Genes showing both causal effects and differential expression
- **Cell-type-specific insights** - Disease mechanisms at cellular resolution
- **Validation results** - Cross-dataset validation of findings
- **Supplementary analyses** - Extended results and validations

## Usage

Run the Crohn's disease analysis scripts in order:

```bash
# Prepare DEG data
Rscript scripts/5-crohns/deg/1-prepare_crohns_deg_pre-processed.R

# Find overlapping genes
Rscript scripts/5-crohns/deg/2-find_deg_mr_overlap_genes_both_colon_and_ti.R

# Prepare comparison data
python scripts/5-crohns/deg/3-prepare_data_for_MR_comparison.py

# Generate visualizations
python scripts/5-crohns/deg/4-plot_UMAP_cell_annotation.py
python scripts/5-crohns/deg/5-plot_UMAP_genes_of_interest_counts.py

# Create figures
Rscript scripts/5-crohns/figures/Figure5-combined_Crohns_figure.R
# ... additional figure scripts as needed
```

**Prerequisites:**
1. Complete Snakemake workflow
2. Run preprocessing and MR analyses
3. Ensure Crohn's disease single-cell data is available

## Key Features

- **Integrated analysis** - Combines genetics, genomics, and single-cell data
- **Cell-type resolution** - Disease mechanisms at cellular level
- **Validation framework** - Multiple orthogonal validation approaches
- **Clinical relevance** - Direct application to human disease

This case study demonstrates the power of integrating genetic and single-cell approaches to understand disease mechanisms and validates the broader analytical framework used throughout the manuscript.