## Purpose: Convert MAGMA gene-level Z statistics to Zscore file for SC-DRS

library(arrow)
library(tidyverse)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

read_parquet(INPUT[[1]]) %>% 
  pivot_wider(id_cols = GENE, names_from = phenotype, values_from = ZSTAT) %>% 
  write_tsv(OUTPUT[[1]])



