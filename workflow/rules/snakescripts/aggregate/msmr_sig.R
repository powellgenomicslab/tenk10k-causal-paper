# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
X <- snakemake@wildcards

gene_universe <- read_lines(INPUT$gene_universe)

q_thresh <- as.numeric(X$q_thresh)
heidi_thresh <- as.numeric(X$heidi_thresh)

df_msmr <- read_parquet(INPUT$msmr) %>% 
  filter(probeID %in% gene_universe,
         b_SMR != 0) %>%
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues, .by = "phenotype") %>% 
  filter(ifelse(is.na(p_HEIDI),
                qval_msmr_pheno < q_thresh,
                qval_msmr_pheno < q_thresh & p_HEIDI >= heidi_thresh))

fwrite(df_msmr, OUTPUT$msmr_sig, sep = "\t", quote = FALSE, row.names = FALSE)