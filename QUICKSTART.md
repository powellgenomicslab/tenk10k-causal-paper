# Quick Start Guide

This guide will help you get up and running with the TenK10K causal inference analysis pipeline quickly.

## Prerequisites

- Git
- Conda or Python 3.8+ with pip
- R 4.0+
- 16+ GB RAM
- 50+ GB free disk space

## 5-Minute Setup

### 1. Clone and Setup Environment
```bash
# Clone repository
git clone https://github.com/powellgenomicslab/tenk10k-causal-paper.git
cd tenk10k-causal-paper

# Setup environment (choose one)
conda env create -f environment.yml && conda activate tenk10k-causal-paper
# OR
pip install -r requirements.txt
```

### 2. Quick Test Run
```bash
# Test with example data (if available)
cd workflow/
snakemake --dry-run --cores 1

# Run preprocessing
Rscript ../scripts/0-preprocess/preprocess_results.R
```

## Full Analysis Pipeline

### Step 1: Data Processing (2-4 hours)
```bash
cd workflow/
snakemake --profile profiles/default --cores 8
```

### Step 2: Generate Results (30-60 minutes)
```bash
cd ..
# Run all analysis scripts
Rscript scripts/0-preprocess/preprocess_results.R
Rscript scripts/1-overview/study_design.R
Rscript scripts/2-mr/mr_results_main.R
Rscript scripts/3-polygenic/scDRS/scdrs_main_supp.R
Rscript scripts/4-drug/drug_enrichment_main.R
```

## Expected Outputs

- **Figures**: Main manuscript figures saved to `figures/`
- **Tables**: Summary statistics and results tables
- **Intermediate data**: Processed datasets in `results/`

## Common Issues

**Memory errors**: Increase available RAM or use smaller datasets
**Missing dependencies**: Check `environment.yml` or `requirements.txt`
**Data not found**: Ensure input data is in expected locations

## Next Steps

1. Read the [full README](README.md) for detailed information
2. Check individual [script documentation](scripts/) for specific analyses
3. See [CONTRIBUTING.md](CONTRIBUTING.md) to contribute improvements

## Getting Help

- Check existing [Issues](https://github.com/powellgenomicslab/tenk10k-causal-paper/issues)
- Create a new issue with the appropriate template
- See individual README files for detailed troubleshooting