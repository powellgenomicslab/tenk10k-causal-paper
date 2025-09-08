# 4-drug: Drug Target Enrichment Analysis

This directory contains scripts for analyzing drug target enrichment in causal genes identified through Mendelian Randomization analyses.

## Overview

This section performs drug target enrichment analyses to identify whether genes with causal effects on disease traits are enriched for known drug targets. This analysis helps bridge genetic discoveries with therapeutic opportunities by connecting causal genes to existing or potential drug interventions.

## Contents

### Scripts
- `drug_enrichment_main.R` - Main drug enrichment analysis script that:
  - Analyzes enrichment of drug targets among causal genes
  - Compares different categories of drug targets
  - Creates visualizations for drug target enrichment results
  - Generates statistical summaries for manuscript

## Dependencies

### R packages
```r
# Standard packages loaded via 0-preprocess/preprocess_results.R
# Additional packages for drug target analysis and visualization
```

## Input Files

The script uses:
- Causal genes identified from MR analyses (from `2-mr/` results)
- Drug target databases and annotations
- Gene annotations and mappings
- Preprocessed results from earlier analysis steps

## Key Analyses

1. **Drug target annotation** - Mapping causal genes to known drug targets
2. **Enrichment testing** - Statistical testing for drug target enrichment
3. **Target category analysis** - Analysis by drug target categories:
   - Approved drug targets
   - Clinical trial targets
   - Investigational targets
4. **Therapeutic area mapping** - Connection to specific therapeutic areas

## Output

Generates:
- **Enrichment statistics** - Statistical significance of drug target enrichment
- **Target gene lists** - Causal genes that are drug targets
- **Visualization plots** - Enrichment plots and target distribution figures
- **Therapeutic insights** - Connections between causal genes and therapeutic opportunities

## Usage

Run the drug enrichment analysis from the repository root:

```bash
Rscript scripts/4-drug/drug_enrichment_main.R
```

**Prerequisites:** 
1. Complete the MR analyses (`2-mr/` scripts)
2. Run the preprocessing script (`0-preprocess/preprocess_results.R`)
3. Ensure drug target databases are available

## Key Features

- **Therapeutic relevance** - Connects genetic discoveries to drug development
- **Multiple target categories** - Comprehensive analysis across target types
- **Statistical validation** - Rigorous statistical testing for enrichment
- **Clinical translation** - Direct relevance to therapeutic development

## Applications

Results from this analysis can inform:
- **Drug repositioning** - Existing drugs for new indications
- **Target prioritization** - Genetic evidence for drug targets
- **Therapeutic development** - Novel targets with causal genetic support
- **Clinical validation** - Genetic evidence supporting drug mechanisms

This analysis provides a crucial link between genetic discoveries and therapeutic applications, helping translate causal genetic insights into potential treatments.