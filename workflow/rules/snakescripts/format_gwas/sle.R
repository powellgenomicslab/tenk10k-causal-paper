library(data.table)
library(tidyverse)

setDTthreads(8)
# setDTthreads(snakemake@threads)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# gwas_file <- "resources/sumstats/gwas/sle.gwas"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"

gwas_df <- fread(gwas_file) %>% 
    filter(chrom %in% 1:22)
gwas_df[, chrom := as.numeric(chrom)]

# Liftover code
source(liftover_script)

coord_df <- gwas_df[, .(chrom, pos, rsid)]
coord_df <- coord_df[, chrom := paste0("chr", chrom)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("chrom", "pos_b38", "rsid")

out_df <- merge(gwas_df, coord_hg38_df, by = c("chrom", "rsid"))

# Get MAF information from TenK10K panel
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2, A1_freq = i.MAF),
       on = c("chrom" = "CHR", "pos_b38" = "POS")]

# format
out_df_format <- out_df %>% 
    filter(!is.na(snp_id)) %>%
    arrange(chrom, pos_b38) %>% 
    mutate(b = ifelse(effect_allele == A1, beta, -beta),
           N = floor(4 / (1/7219 + 1/15991))) %>% 
    select(SNP = snp_id, A1, A2, freq = A1_freq, b, se,  p, N)


# Write results
out_file <- snakemake@output[[1]]
# out_file <- "resources/pipeline_ma/sle.ma"
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)
