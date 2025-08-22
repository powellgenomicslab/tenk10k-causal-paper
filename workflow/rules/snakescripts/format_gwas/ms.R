# summary stats downloaded from meta-analysis (finngen r12 - mvp - ukbb)
# Manuscript: https://www.science.org/doi/10.1126/science.aav7188
# according to the supplementary table from the published paper,
# A1 seems to be the effect allele
# N case = 47,429
# N control = 68,374
library(data.table)
library(tidyverse)

setDTthreads(snakemake@threads)

# setDTthreads(4)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# gwas_file <- "resources/sumstats/gwas/ms.gwas"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"

gwas_df <- fread(gwas_file) %>% 
    filter(!is.na(P))

# liftover
source(liftover_script)

coord_df <- gwas_df[, .(CHR, BP, SNP)] 
coord_df <- coord_df[, CHR := paste0("chr", CHR)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("CHR", "pos_b38", "SNP")

# Merge
out_df <- merge(gwas_df, coord_hg38_df, by = c("CHR", "SNP"))

# get allele frequency from tenk10k since there is no default
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

setnames(out_df, c("A1", "A2"), c("EA", "OA"))

out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2, A1_freq = i.MAF),
       on = c("CHR", "pos_b38" = "POS")]

# format
out_df_format <- out_df %>% 
    filter(!is.na(A1_freq)) %>%
    arrange(CHR, pos_b38) %>% 
    mutate(abs_z = qnorm(1 - P/2),
           b = ifelse(EA == A1, log(OR), -log(OR)),
           se = abs(b) / abs_z,
           N = floor(4 / (1/47429 + 1/68374))) %>% 
    select(SNP = snp_id, A1, A2, freq = A1_freq, b, se,  p = P, N)

# Write results
out_file <- snakemake@output[[1]]
# out_file <- "resources/pipeline_ma/ms.ma"
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)