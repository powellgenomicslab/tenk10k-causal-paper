# Script for processing summary statistics from the IPF study by GBMI
# Column name	Description
# #CHR	chromosome
# POS	genome coordinate (hg38)
# REF	reference allele
# ALT	alternative allele
# rsid	rs id
# all_meta_AF	frequency of ALT in the meta-analysis
# inv_var_meta_beta	effect size of ALT in the meta-analysis
# inv_var_meta_sebeta	standard error for effect size
# inv_var_meta_p	association p value in the meta-analysis
# inv_var_het_p	heterogeneity (Cochran’s Q test) p value
# direction	direction of effect sizes in all data sets in the meta-analysis
# N_case	number of cases
# N_ctrl	number of controls
# n_bbk	number of biobanks in the meta-analysis (minimum is 2)
# is_strand_flip	flag to indicate whether this variant had a strand flip in any data set
# is_diff_AF_gnomAD	flag to indicate whether this variant failed the AF QC (compared to gnomAD) in any data set

library(data.table)

setDTthreads(snakemake@threads)

# Read in summary statistics
gwas_df <- fread(snakemake@input[[1]])

colnames(gwas_df)[1] <- "CHR"
gwas_df <- gwas_df[CHR %in% 1:22]
gwas_df <- gwas_df[,CHR := as.numeric(CHR)]

# Set SNP ids
gwas_df[, A1 := fifelse(all_meta_AF < 0.5, ALT, REF)]
gwas_df[, A2 := fifelse(all_meta_AF < 0.5, REF, ALT)]
gwas_df[, freq := fifelse(all_meta_AF < 0.5, all_meta_AF, 1 - all_meta_AF)]
gwas_df[, beta := fifelse(all_meta_AF < 0.5, inv_var_meta_beta, -1 * inv_var_meta_beta)]
gwas_df[, SNP := paste0(CHR, ":", POS, ":", A2, ":", A1)]

# Select columns to output
old_columns <-  c("SNP", "A1", "A2", "freq", "beta", "inv_var_meta_sebeta", "inv_var_meta_p", "n_dataset")
new_columns <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
output_df <- gwas_df[, ..old_columns]
colnames(output_df) <- new_columns

# Write results    
fwrite(output_df, snakemake@output[[1]], sep = "\t", na = "NA", quote = FALSE)
