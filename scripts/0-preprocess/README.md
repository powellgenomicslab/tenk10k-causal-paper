# 0-preprocess: Data Preprocessing

This directory contains scripts and data for preprocessing GWAS summary statistics and eQTL data for downstream causal inference analyses.

## Overview

The preprocessing step is crucial for ensuring data quality and consistency across different datasets used in the causal inference pipeline. This includes formatting GWAS summary statistics, preparing eQTL data, and creating unified gene and trait mappings.

## Contents

### Scripts
- `preprocess_results.R` - Main preprocessing script that:
  - Loads cell type and trait mappings from metadata
  - Processes MAGMA gene-level association results
  - Filters and harmonizes multi-SMR (mSMR) results from TenK10K and eQTLGen
  - Creates a unified gene universe for downstream analyses
  - Applies quality control filters (removes missing p-values, zero effect sizes)

### Subdirectories
- `gwas/` - Contains preprocessing scripts and metadata for GWAS summary statistics
- `tenk10k-eqtl/` - Contains preprocessing scripts for TenK10K Phase 1 single-cell eQTL data

## Dependencies

### R packages
```r
library(data.table)    # Fast data manipulation
library(tidyverse)     # Data wrangling and visualization
library(arrow)         # Parquet file handling
library(fs)            # File system operations
library(qvalue)        # Multiple testing correction
```

## Input Files

The script expects the following input files in specific locations:
- `metadata/cell_map.tsv` - Cell type mapping information
- `metadata/trait_map.tsv` - Trait mapping and inclusion criteria
- `metadata/trait_category.tsv` - Trait category definitions
- `resources/misc/gencode.v44.gene_type.tsv` - Gene annotations
- `results/aggregate/tenk10k_phase1.magma.gz.parquet` - MAGMA results
- `results/aggregate/tenk10k_phase1.msmr.gz.parquet` - TenK10K mSMR results
- `results/aggregate/eqtlgen2020.msmr.gz.parquet` - eQTLGen mSMR results

## Output

The preprocessing generates filtered and harmonized datasets ready for downstream causal inference analyses, including:
- Quality-controlled gene-level association results
- Unified gene universe across different datasets
- Standardized trait and cell type mappings

## Usage

Run the preprocessing script from the repository root:

```bash
Rscript scripts/0-preprocess/preprocess_results.R
```

**Note:** Ensure that the Snakemake workflow has been executed first to generate the required input files in `results/aggregate/`.

## Quality Control

The script applies several QC filters:
- Removes genes with missing SMR p-values
- Excludes associations with zero GWAS or SMR effect sizes
- Creates intersection of genes available in both MAGMA and TenK10K analyses
- Applies trait inclusion criteria from metadata