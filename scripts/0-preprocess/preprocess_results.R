# preprocess main results

library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(qvalue)

df_cell_map <- fread("metadata/cell_map.tsv")
df_trait_map_all <- fread("metadata/trait_map.tsv")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

setDT(df_trait_map)

df_trait_cat <- fread("metadata/trait_category.tsv")
cat_order <- df_trait_cat$cat_order

df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")

df_msmr_tenk10k <- read_parquet("results/aggregate/tenk10k_phase1.msmr.gz.parquet") %>%
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0)

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
gene_universe <- intersect(df_magma_all$GENE, df_msmr_tenk10k$probeID)

df_msmr_eqtlgen <- read_parquet("results/aggregate/eqtlgen2020.msmr.gz.parquet") %>%
  filter(
    !is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0,
    probeID %in% gene_universe
  ) %>%
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues, .by = "phenotype") %>%
  mutate(sig = ifelse(is.na(p_HEIDI), qval_msmr_pheno < 0.05,
    qval_msmr_pheno < 0.05 & p_HEIDI >= 0.01
  )) %>%
  setDT()

phenotypes <- df_trait_map %>%
  pull(trait_id)

df_magma <- df_magma_all %>%
  filter(GENE %in% gene_universe, phenotype %in% phenotypes) %>%
  mutate(qval = qvalue(P)$qvalues) %>%
  filter(qval < 0.05) %>%
  setDT()

df_msmr_tenk10k <- df_msmr_tenk10k %>%
  filter(probeID %in% gene_universe) %>%
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues, .by = "phenotype") %>%
  mutate(sig = ifelse(is.na(p_HEIDI), qval_msmr_pheno < 0.05,
    qval_msmr_pheno < 0.05 & p_HEIDI >= 0.01
  )) %>%
  left_join(df_cell_map %>% select(biosample = wg2_scpred_prediction, cell_type, major_cell_type)) %>%
  inner_join(df_trait_map %>%
    select(
      phenotype = trait_id, pheno_label = label,
      pheno_cat = cat_rev, supercategory
    )) %>%
  group_by(phenotype) %>%
  mutate(
    cell_type = factor(cell_type, df_cell_map$cell_type),
    major_cell_type = factor(major_cell_type, unique(df_cell_map$major_cell_type)),
    pheno_cat = factor(pheno_cat, cat_order)
  ) %>%
  setDT()

# annotate gene and phenotype
df_msmr_tenk10k[df_gene_annot, gene_type := i.gene_type, on = c("probeID" = "ensembl_gene_id")]
df_msmr_tenk10k[, magma_gene := FALSE]
df_msmr_tenk10k[df_magma, magma_gene := TRUE, on = c("probeID" = "GENE", "phenotype")]
df_msmr_tenk10k[, eqtlgen_mr := FALSE]
df_msmr_tenk10k[df_msmr_eqtlgen, eqtlgen_mr := i.sig, on = c("probeID", "phenotype")]
# filter to only significant gene
df_msmr <- df_msmr_tenk10k[sig == TRUE]

pheno_order <- df_msmr[, .N, by = pheno_label][order(-N), pheno_label]
df_msmr[, pheno_label := factor(pheno_label, pheno_order)]
