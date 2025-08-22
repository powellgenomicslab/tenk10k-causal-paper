##########################################################
#
# Processing summary statistics for migraine
# DOI: https://doi.org/10.1038/s41588-021-00990-0
# Reference: hg38
# Reference allele is the effect allele
# Other allele is the reference allele
# 
#########################################################

library(data.table)
library(tidyverse)

setDTthreads(8)
# setDTthreads(snakemake@threads)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# gwas_file <- "resources/sumstats/gwas/migraine.gwas"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"

gwas_df <- fread(gwas_file) %>% 
    filter(chromosome %in% 1:22)
gwas_df[, chromosome := as.numeric(chromosome)]

# Liftover code
source(liftover_script)

coord_df <- gwas_df[, .(chromosome, position, rs_number)] 
coord_df <- coord_df[, chromosome := paste0("chr", chromosome)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("chromosome", "pos_b38", "rs_number")

# Merge
out_df <- merge(gwas_df, coord_hg38_df, by = c("chromosome", "rs_number"))

# MAF check
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2, A1_freq = i.MAF),
       on = c("chromosome" = "CHR", "pos_b38" = "POS")]

# Select columns to output
out_df_format <- out_df %>% 
    filter(!is.na(snp_id)) %>%
    arrange(chromosome, pos_b38) %>% 
    mutate(b = ifelse(reference_allele == A1, beta, -beta),
           freq = ifelse(reference_allele == A1, eaf, 1-eaf),
           maf_calc = pmin(eaf, 1-eaf)) %>% 
    # remove variants where MAF diff >= 0.2
    filter(abs(A1_freq - maf_calc) < 0.2)

# infer indels based on allele frequency
out_df_format[reference_allele %in% c("D", "I"),
              `:=`(b = ifelse(eaf < 0.5, beta, -beta),
                   freq = ifelse(eaf < 0.5, eaf, 1-eaf))]

out_df_format <- out_df_format %>% 
    select(SNP = snp_id, A1, A2, freq, b, se, p = p.value, N = Neff) 

# Write results
out_file <- snakemake@output[[1]]
# out_file <- "resources/pipeline_ma/migraine.ma"
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)
