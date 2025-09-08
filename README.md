# TenK10K Phase1 causal inference manuscript

![Static Badge](https://img.shields.io/badge/version-preprint-red)

> [!NOTE]
> This repository is still under active development and will be updated from time to time

## Study design
![](figures/biorender/study_design.png)

## Overview

This repository contains code, data, and workflows for the TenK10K causal inference manuscript. The project is organized into several directories, each corresponding to a major analysis step or component.

## Directory Structure

- **`figures/`**  
  Contains figures generated for the manuscript, including study design diagrams.

- **`metadata/`**  
  Contains metadata files for traits, trait categories, and cell types from TenK10K Phase 1 sc-eQTL analysis.

- **`scripts/`**  
  Main analysis scripts, organized by analysis section:
  - **`0-preprocess/`** - Data preprocessing, GWAS and eQTL data preparation
  - **`1-overview/`** - Study overview and summary statistics  
  - **`2-mr/`** - Mendelian Randomization analyses and comparisons
  - **`3-polygenic/`** - Polygenic enrichment analyses (scDRS, scDeepID)
  - **`4-drug/`** - Drug target enrichment and therapeutic relevance
  - **`5-crohns/`** - Crohn's disease case study with matched single-cell data
  - **`util/`** - Utility functions and helper scripts

- **`workflow/`**  
  Snakemake pipeline for reproducible data processing and analysis. The workflow handles data formatting, quality control, statistical analyses, and intermediate file generation.

## Key Features

- **Reproducible analysis pipeline** using Snakemake workflow management
- **Multi-scale integration** of GWAS, eQTL, and single-cell data
- **Cell-type-specific causal inference** using TenK10K Phase 1 data
- **Comprehensive validation** with external datasets (eQTLGen, MAGMA)
- **Therapeutic translation** through drug target enrichment analysis
- **Disease case study** demonstrating framework application

## Methods Overview

This repository implements several key analytical approaches:

1. **Summary-based Mendelian Randomization (SMR)** - Tests for causal effects of gene expression on complex traits
2. **Multi-SMR (mSMR)** - Joint analysis across multiple cell types
3. **scDRS** - Single-cell disease relevance scoring for polygenic enrichment
4. **Drug target enrichment** - Assessment of therapeutic relevance for causal genes
5. **Integration analysis** - Cross-validation with bulk tissue and gene-level results

## Installation and Setup

### 1. Clone the repository
```bash
git clone https://github.com/powellgenomicslab/tenk10k-causal-paper.git
cd tenk10k-causal-paper
```

### 2. Set up the environment

#### Option A: Using Conda (Recommended)
```bash
# Create environment from file
conda env create -f environment.yml

# Activate the environment
conda activate tenk10k-causal-paper
```

#### Option B: Using pip for Python dependencies only
```bash
# Install Python dependencies
pip install -r requirements.txt

# Note: R dependencies need to be installed separately
```

#### Option C: Manual installation
Install the following dependencies manually:

**R packages:**
```r
install.packages(c("tidyverse", "data.table", "arrow", "patchwork", 
                   "scales", "paletteer", "ragg", "readxl", "writexl",
                   "broom", "fs", "glue", "ggrepel", "ggridges", 
                   "formattable", "ggnewscale", "geomtextpath", "ggforce"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("qvalue", "rtracklayer"))
```

**Python packages:**
```bash
pip install pandas numpy matplotlib seaborn scanpy torch scipy scikit-learn snakemake scdrs gget
```

**External tools:**
- MAGMA (for gene-level association analysis)
- PLINK (for genetic data processing)
- CrossMap (for genomic coordinate conversion)

### 3. Data Requirements

Before running analyses, ensure you have access to:
- GWAS summary statistics (formatted as specified in `workflow/README.md`)
- TenK10K Phase 1 single-cell eQTL results
- Reference genetic data (LD reference panels, gene annotations)

## Usage

### Quick Start

1. **Run the Snakemake workflow** (generates intermediate results):
   ```bash
   cd workflow/
   snakemake --profile profiles/default --cores 8
   ```

2. **Run analysis scripts** (generates figures and tables):
   ```bash
   # Preprocessing
   Rscript scripts/0-preprocess/preprocess_results.R
   
   # Overview analysis
   Rscript scripts/1-overview/study_design.R
   
   # Mendelian Randomization
   Rscript scripts/2-mr/mr_results_main.R
   
   # Polygenic analyses
   Rscript scripts/3-polygenic/scDRS/scdrs_main_supp.R
   
   # Drug target analysis
   Rscript scripts/4-drug/drug_enrichment_main.R
   
   # Crohn's disease case study
   Rscript scripts/5-crohns/figures/Figure5-combined_Crohns_figure.R
   ```

### Detailed Instructions

For detailed instructions on each analysis step, see the README files in each subdirectory:
- [`workflow/README.md`](workflow/README.md) - Snakemake pipeline setup and execution
- [`scripts/0-preprocess/README.md`](scripts/0-preprocess/README.md) - Data preprocessing
- [`scripts/1-overview/README.md`](scripts/1-overview/README.md) - Study overview and design
- [`scripts/2-mr/README.md`](scripts/2-mr/README.md) - Mendelian Randomization analysis
- [`scripts/3-polygenic/README.md`](scripts/3-polygenic/README.md) - Polygenic enrichment analysis
- [`scripts/4-drug/README.md`](scripts/4-drug/README.md) - Drug target enrichment
- [`scripts/5-crohns/README.md`](scripts/5-crohns/README.md) - Crohn's disease case study

## Execution Order

The analyses should be run in the following order:

1. **Snakemake workflow** (`workflow/`) - Generates intermediate data files
2. **Preprocessing** (`scripts/0-preprocess/`) - Prepares unified datasets  
3. **Analysis scripts** (`scripts/1-overview/` through `scripts/5-crohns/`) - Can be run in parallel or sequentially

## Expected Runtime

- **Snakemake workflow**: 2-6 hours (depending on data size and computational resources)
- **Analysis scripts**: 30 minutes - 2 hours total
- **Full pipeline**: 3-8 hours total

## Hardware Requirements

- **Memory**: 16+ GB RAM recommended (32+ GB for large datasets)
- **Storage**: 50+ GB free space for intermediate files
- **CPU**: 4+ cores recommended for parallel processing

## Troubleshooting

- Check individual README files for script-specific troubleshooting
- Ensure all dependencies are properly installed
- Verify data file paths and formats
- Check the [Issues](https://github.com/powellgenomicslab/tenk10k-causal-paper/issues) page for known problems

## Data Availability

- **Single-cell eQTL data**: TenK10K Phase 1 results (access information to be updated)
- **GWAS summary statistics**: Links to original data sources provided in workflow documentation
- **Processed results**: Key analysis outputs will be made available upon publication

## Contributing

We welcome contributions to improve the analysis pipeline or documentation. Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

For bug reports or feature requests, please open an [issue](https://github.com/powellgenomicslab/tenk10k-causal-paper/issues).

## Citation

If you use this code or data in your research, please see [`CITATION.md`](CITATION.md) for citation information.

## License

This project is licensed under the MIT License - see the [`LICENSE`](LICENSE) file for details.

## Contact

For questions about this analysis or collaboration opportunities:
- Open an issue on GitHub: https://github.com/powellgenomicslab/tenk10k-causal-paper/issues
- Contact the Powell Genomics Lab: [Contact information to be updated]

## Acknowledgments

- TenK10K Consortium members and contributors
- Data providers and consortia (eQTLGen, GWAS consortia)
- Software developers (Snakemake, R/Bioconductor, Python scientific computing ecosystem)
- Computational resources and support teams
