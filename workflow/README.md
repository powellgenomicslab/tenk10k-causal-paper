# Snakemake Workflow for Tenk10k Causal Paper

## Overview
This Snakemake workflow processes and analyzes data for the Tenk10k causal inference paper, ensuring reproducibility and scalability.

## Features
- Data preprocessing
- Liftover of GWAS summary statistics
- Gene set enrichment analysis (gget Enrichr)
- SMR analysis
- Integration with Open Targets Platform
- MAGMA-based gene-level association analysis
- scDRS computation

## Requirements
- Python 3.x
- Snakemake
- R (v4.0+) with required libraries
- Conda (optional)
- External tools:
  - MAGMA
  - PLINK
  - CrossMap

## Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/tenk10k-workflow.git
   cd tenk10k-workflow
   ```
2. Install dependencies:
   - Using Conda:
     ```bash
     conda env create -f environment.yaml
     conda activate tenk10k
     ```
   - Or install Python and R packages manually.
3. Configure paths in `config/finngen_meta_path.yaml`.

## Usage
Run the workflow:
```bash
snakemake --profile profiles/default
```

### Example Commands

- Dry-run the workflow:
```bash
snakemake --profile profiles/default --dry-run
```

- Run with specific targets:
```bash
snakemake --profile profiles/default results/enrichment/{study}.tsv
```

- Use multiple cores:
```bash
snakemake --profile profiles/default --cores 8
```

## Workflow Steps

### 1. GWAS Formatting
- **Scripts:** `rules/snakescripts/format_gwas/*.R`
- **Purpose:** Standardize GWAS summary statistics, liftover coordinates, harmonize alleles, filter variants.
- **Output:** Formatted files in `resources/pipeline_ma/`

### 2. eQTL Preparation
- **Scripts:** `rules/snakescripts/prep_besd_chr/*.sh`, `prep_smr_input/*.R`
- **Purpose:** Prepare BESD files, probe lists, and p-value thresholds for SMR.
- **Output:** `resources/besd/`, `resources/smr/`

### 3. SMR Analysis
- **Scripts:** `rules/snakescripts/run_smr.sh`, `smr_locus.sh`
- **Purpose:** Run SMR per chromosome/cell type/phenotype, extract locus information.
- **Output:** `results/smr/`, `results/smr_locus/`

### 4. MAGMA Gene-Level Association
- **Scripts:** `rules/snakescripts/aggregate/magma.R`, `magma_format_output.R`
- **Purpose:** Aggregate MAGMA results, annotate genes, FDR correction.
- **Output:** `results/magma/`

### 5. Gene Set Enrichment (Enrichr)
- **Scripts:** `rules/snakescripts/enrichment/gget_enrichr.py`, `gget_enrichr_pheno.py`
- **Purpose:** Perform gene set enrichment using Enrichr via gget.
- **Output:** `results/enrichment/`

### 6. scDRS Computation
- **Scripts:** `rules/snakescripts/scdrs/*`
- **Purpose:** Prepare covariates, regress out confounders, compute scDRS scores.
- **Output:** `resources/scdrs/`, `results/scdrs/`

### 7. Genetic Correlation
- **Scripts:** `rules/snakescripts/ldak/*`
- **Purpose:** Calculate genetic correlations between traits using LDAK.
- **Output:** `results/gen_cor/`
