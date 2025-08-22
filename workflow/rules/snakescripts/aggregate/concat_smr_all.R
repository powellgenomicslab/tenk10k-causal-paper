# Concatenate SMR results from a study

library(tidyverse)
library(arrow)
library(glue)
library(fs)
library(qvalue)

files <- snakemake@input
PARAMS <- snakemake@params

# PARAMS <- list(
#     fdr_msmr = 0.05
# ) 

read_data <- function(x){
    file_parts <- path_split(x)[[1]]
    biosample <- file_parts[4]
    phenotype <- file_parts[5]

    data <- read_tsv_arrow(x) %>%
        mutate(biosample = biosample,
               phenotype = phenotype,
               .before = probeID)
}

calc_q <- function(p, ...) {
    possibly(qvalue, qvalue_truncp(p))(p, ...) %>% .$qvalues
}

df_all <- map_df(files, read_data)  %>% 
    mutate(qval_msmr_biosample_pheno = calc_q(p_SMR_multi), .by = c("biosample", "phenotype")) %>% 
    mutate(qval_msmr_biosample = calc_q(p_SMR_multi), .by = c("biosample")) %>% 
    mutate(qval_msmr_pheno = calc_q(p_SMR_multi), .by = c("phenotype")) %>% 
    mutate(qval_msmr = calc_q(p_SMR_multi))

# write results
write_parquet(
    df_all,
    snakemake@output[[1]],
    compression = "gzip"
)

# qvalue adjustment

# write results

# df_all <- df_all %>% 
#     group_by(biosample) %>% 
#     mutate(qval_msmr = qvalue(p_SMR_multi)$qvalues,
#            qval_singlesmr = qvalue(p_SMR)$qvalues)