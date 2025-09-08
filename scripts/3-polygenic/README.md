# 3-polygenic: Polygenic Enrichment Analysis

This directory contains scripts for polygenic enrichment analyses using single-cell Disease Relevance Score (scDRS) and scDeepID methods to identify cell-type-specific disease associations.

## Overview

This section performs polygenic enrichment analyses to identify which cell types show enrichment for disease-associated genetic variants. The analysis integrates two complementary approaches: scDRS for disease relevance scoring and scDeepID for cell type identification, providing insights into cell-type-specific disease mechanisms.

## Contents

### Subdirectories

#### `scDRS/`
Contains scripts for single-cell Disease Relevance Score analysis:
- `scdrs_main_supp.R` - Main scDRS analysis and figures for manuscript
- `scdrs_supp_sampling.R` - Supplementary sampling analyses for scDRS
- `README.md` - Detailed documentation for scDRS analyses

#### `scDeepID/` 
Contains scripts for scDeepID integration analysis:
- `scDRS_scDeepID_BCELL.py` - Integration analysis focusing on B cell populations
- `README.md` - Documentation for scDRS-scDeepID integration

## Methods

### scDRS (single-cell Disease Relevance Score)
- Calculates disease relevance scores for individual cells
- Identifies cell types enriched for disease-associated variants
- Provides statistical framework for cell-type-specific disease associations
- Integrates GWAS summary statistics with single-cell expression data

### scDeepID Integration
- Combines scDRS results with scDeepID cell type predictions
- Focuses on specific cell populations (e.g., B cells)
- Validates disease associations across different cell type annotation methods
- Provides orthogonal validation of cell-type-specific effects

## Dependencies

### R packages (scDRS)
```r
# Standard packages loaded via 0-preprocess/preprocess_results.R
# Additional visualization and analysis packages as needed
```

### Python packages (scDeepID)
```python
# Single-cell analysis packages
# Integration with scDRS outputs
# Statistical analysis and visualization
```

## Input Files

The scripts use:
- scDRS score results from the Snakemake workflow
- scDeepID cell type predictions
- GWAS summary statistics
- Single-cell expression data from TenK10K Phase 1

## Key Analyses

1. **Disease relevance scoring** - Cell-level disease association scores
2. **Cell type enrichment** - Statistical testing for cell-type-specific enrichment
3. **Cross-method validation** - Comparison between scDRS and scDeepID results
4. **Population-specific analysis** - Focused analysis on specific cell populations

## Output

Generates:
- **Cell-type enrichment results** - Statistical significance of disease associations
- **Visualization plots** - Heatmaps, scatter plots, and enrichment visualizations
- **Integration analyses** - Combined scDRS-scDeepID results
- **Population-specific insights** - Detailed analysis of key cell types

## Usage

Run the polygenic analysis scripts from the repository root:

```bash
# scDRS main analysis
Rscript scripts/3-polygenic/scDRS/scdrs_main_supp.R

# scDRS sampling analysis
Rscript scripts/3-polygenic/scDRS/scdrs_supp_sampling.R

# scDeepID integration (B cell focus)
python scripts/3-polygenic/scDeepID/scDRS_scDeepID_BCELL.py
```

**Prerequisites:** 
1. Complete the Snakemake workflow to generate scDRS scores
2. Run the preprocessing script (`0-preprocess/preprocess_results.R`)
3. Ensure scDeepID results are available for integration analyses

## Key Features

- **Multi-method approach** - Combines scDRS and scDeepID for robust results
- **Cell-type specificity** - Identifies disease-relevant cell populations
- **Statistical rigor** - Comprehensive statistical testing and validation
- **Integration focus** - Cross-validates results across methods

See individual README files in `scDRS/` and `scDeepID/` subdirectories for detailed method-specific documentation.