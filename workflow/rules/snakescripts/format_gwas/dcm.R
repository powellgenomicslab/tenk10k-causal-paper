# Script for processing summary statistics from the HF study by Albert Henry
# Column Description
#key Unique variant ID in <chr>:<pos_b37>:<A1>_<A2> format
# rsID dbSNP rsID
# chr chromosome
# pos_b37 base pair coordinate (genome build b37 / hg19)
# A1 effect allele (first allele in alphabetical order)
# A2 non-effect allele (second allele in alphabetical order)
# A1_beta
# log odds ratio estimated from fixed-effect meta-analysis implemented in METAL per one
# allele increase of A1
# A1_freq effect allele frequency
# se standard error of A1_beta
# pval P-value for association with phenotype
# logP -log10 P-value for association with phenotype
# N_case variant-specific number of cases contributed to the meta-analysis
# N_total variant-specific total number of sample (case + control) contributed to the meta-analysis
# isq_het I squared statistic for heterogeneity estimated using METAL
# p_het P-value for heterogeneity estimated using METAL

# To do - liftover to hg38

library(data.table)
library(readxl)
library(tidyverse)

setDTthreads(snakemake@threads)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# gwas_file <- "resources/sumstats/gwas/dcm.gwas"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"

source(liftover_script)

# Read in summary statistics
gwas_df <- fread(gwas_file)

# Liftover

coord_df <- gwas_df[, .(chr, pos_b37, `#key`)] 
coord_df <- coord_df[, chr := paste0("chr", chr)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("chr", "pos_b38", "#key")

# Merge
gwas_df <- merge(gwas_df, coord_hg38_df, by = c("chr", "#key"))

# Select columns to output
gwas_df <- gwas_df[chr %in% 1:22]

# Flip results by effect allele being defined as the minor allele frequency being less than 0.5
setnames(gwas_df, c("A1", "A2"), c("EA", "OA"))
gwas_df[, A1 := fifelse(A1_freq < 0.5, EA, OA)]
gwas_df[, A2 := fifelse(A1_freq < 0.5, OA, EA)]
gwas_df[, freq := fifelse(A1_freq < 0.5, A1_freq, 1 - A1_freq)]
gwas_df[, beta := fifelse(A1_freq < 0.5, A1_beta, -1 * A1_beta)]
gwas_df[, SNP := paste0(chr, ":", pos_b38, ":", A2, ":", A1)]

# Select columns to output
old_columns <-  c("SNP", "A1", "A2", "freq", "beta", "se", "pval", "N_total")
new_columns <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
output_df <- gwas_df[, ..old_columns]
colnames(output_df) <- new_columns

# Write results  
output_file <- snakemake@output[[1]]
# output_file <- "resources/pipeline_ma/dcm.ma"  
fwrite(output_df, output_file, sep = "\t", na = "NA", quote = FALSE)