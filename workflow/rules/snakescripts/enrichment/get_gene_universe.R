## Purpose: Generate gene universe by intersecting MAGMA and MSMR gene sets
## Input: MAGMA and MSMR gene-level .parquet files
## Output: Gene universe text file

library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
df_magma <- read_parquet(INPUT$magma)
df_msmr <- read_parquet(INPUT$msmr)

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
gene_universe <- intersect(df_magma$GENE, df_msmr$probeID)

writeLines(gene_universe, OUTPUT$gene_universe)
