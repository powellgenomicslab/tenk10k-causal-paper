# R scripts to format RA GWAS summary statistics
# GWAS: Ishigaki et al (2022) Nature Genetics
# File: - EUR_all_auto-10-2021.txt.gz: EUR meta-analysis, autosomal (25 cohorts: 22,350 cases and 74,823 controls)
# [File format (header)]
# - SNP: chr_pos_ref_alt (hg19)
# - Beta: effect size estimates (the effect allele is the alternative allele)
# - SE: S.E. of the effect size estimate
# - Pval: P value

library(data.table)
library(tidyverse)

setDTthreads(snakemake@threads)
# setDTthreads(4)

gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]
liftover_script <- snakemake@input[["liftover_script"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"
# gwas_file <- "resources/sumstats/gwas/ra.gwas"

source(liftover_script)
df <- fread(gwas_file) %>% 
    separate(SNP, into = c("chr", "pos_b37", "ref", "alt"), sep = "_", remove = FALSE)

setDT(df)

# liftover
coord_df <- df[, .(chr, pos_b37, SNP)] 
coord_df[, chr := paste0("chr", chr)]

setnames(coord_df, c("seqnames", "start", "snpid"))

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("chr", "pos_b38", "SNP")

# Merge back
df[, chr := as.numeric(chr)]
out_df <- merge(df, coord_hg38_df, by = c("chr", "SNP")) %>% 
    filter(chr %in% 1:22)  %>% 
    arrange(chr, pos_b38)

# get allele frequency from tenk10k since there is no default
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2, A1_freq = i.MAF),
       on = c("chr" = "CHR", "pos_b38" = "POS")]

# Select columns to output
out_df_format <- out_df %>% 
    filter(!is.na(A1_freq)) %>%
    mutate(b = ifelse(alt == A1, Beta, -Beta),
            # N effective based on: - EUR_all_auto-10-2021.txt.gz: EUR meta-analysis, autosomal (25 cohorts: 22,350 cases and 74,823 controls)
           N = floor(4 / (1/22350 + 1/74823))) %>% 
    select(SNP = snp_id, A1, A2, freq = A1_freq, b, se = SE, p = Pval, N) 

# Write results
out_file <- snakemake@output[[1]]
# out_file <- "resources/pipeline_ma/ra.ma"
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)
