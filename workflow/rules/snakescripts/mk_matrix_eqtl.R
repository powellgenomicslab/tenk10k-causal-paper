# snakemake script to create a matrix eQTL file for SMR

library(data.table)
library(tidyverse)
library(fs)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params
THREADS <- snakemake@threads

setDTthreads(THREADS)

args = commandArgs(trailingOnly=TRUE)
cell_type <- args[1]
#cell_type <- "CD4_TCM"
#chr <- "21"

library(data.table)

# Directories
proj_dir <- "/g/data1a/ei56/as8574/analysis/TenK10K_SMR/"
source_data_dir <- paste0(proj_dir, "source_files/")
saige_eqtl_dir <- paste0(source_data_dir, "SAIGE-QTL/2024-12_Freeze/", cell_type, "/")
matrix_eqtl_dir <- paste0(source_data_dir, "MatrixEQTL/2024-12_Freeze/", cell_type, "/")
smr_input_dir <- paste0(proj_dir, "inputs/")
#plink_dir <- paste0("/g/data/fy54/genotypes/common_variants/2024-12_Freeze/plink/")

# Input files
saigeqtl_filepath1 <- paste0(saige_eqtl_dir, cell_type, "_common_all_cis_raw_pvalues.tsv")
saigeqtl_filepath2 <- paste0(saige_eqtl_dir, cell_type, "_all_cis_cv_gene_level_results.tsv")
input_df1 <- fread(INPUT$common_raw)
input_df2 <- fread(INPUT$common_gene)

# Remove results with maf < 1%
input_df1 <- input_df1[pmin(AF_Allele2, 1-AF_Allele2) >= 0.01]

# Adjust T-stat for MatrixEQTL
# Calculate sig results

# qvalue_obj <- qvalue(input_df2$ACAT_p)
# input_df2 <- input_df2[, `:=` (ACAT_qvalue = qvalue_obj$qvalues, ACAT_localFDR = qvalue_obj$lfdr)]
# input_df2 <- input_df2[ACAT_qvalue < 0.05]

# # Get CD4_Naive nominal p-value
# threshold <- max(input_df2$top_pval)

# # Write out results
# threshold_out <- data.frame(cell_id = cell_type, threshold = threshold)
# threshold_filename <- paste0(smr_input_dir, "thresholds/", cell_type, "_Thresholds.csv")
# write.table(threshold_out, threshold_filename, col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)

# for (x in seq(1,22)){
#     subset_df <- input_df1[CHR == x]
#     chr_qobj <- qvalue(subset_df$`p.value`, pi0 = 1)
#     subset_df <- subset_df[, FDR := chr_qobj$qvalue]
#     subset_df[, Tstat := BETA / SE]

#     # Crosscheck with PLINK file
#     # Variants match, should be ok
#     # bim_df <- fread(paste0(plink_dir, "chr", x, "_common_variants.bim"), header)

#     # Need to reformat the SNPID 
#     matrix_df <- subset_df[ , .(MarkerID, gene, BETA, Tstat, `p.value`, FDR)]

#     colnames(matrix_df) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
#     matrix_eqtl_filename <-paste0(matrix_eqtl_dir, cell_type, "_Chr", x, "_MatrixEQTL.tsv")
#     fwrite(matrix_df, matrix_eqtl_filename, sep = "\t", na = "NA", quote = FALSE)
# }