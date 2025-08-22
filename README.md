# TenK10K Phase1 causal inference manuscript

![Static Badge](https://img.shields.io/badge/version-preprint-red)


## Study design
![](figures/biorender/study_design.png)

## Overview

This repository contains code, data, and workflows for the TenK10K causal inference manuscript. The project is organized into several directories, each corresponding to a major analysis step or component.

## Directory Structure

- `figures/`  
  Contains figures generated for the manuscript.

- `metadata/`  
  Contains metadata for traits, trait categories, and cell types from TenK10K phase 1 sc-eQTL analysis.

- `scripts/`  
  Main analysis scripts, organized by section:
  - `0-preprocess/`  
    Preprocessing scripts, including GWAS and eQTL data preparation.
  - `1-overview/`  
    Overview and summary scripts.
  - `2-mr/`  
    Mendelian Randomization analyses.
  - `3-polygenic/`  
    Polygenic enrichment analyses, including scDeepID and scDRS.
  - `4-drug/`  
    Drug target support analyses.
  - `5-crohns/`  
    Crohn's disease-specific analyses.

- `workflow/`  
  Contains the Snakemake pipeline for reproducible analysis. The workflow orchestrates data processing and analysis steps, producing outputs that are further processed by scripts.

## Getting Started

1. Clone the repository:
   ```sh
   git clone https://github.com/powellgenomicslab/tenk10k-causal-paper.git
   cd tenk10k-causal-paper
   ```

2. Run the Snakemake workflow in `workflow/` to generate analysis outputs.

3. Use the scripts in `scripts/` to produce figures and tables for the manuscript.

4. See individual README.md files in each subdirectory for step-specific instructions.
