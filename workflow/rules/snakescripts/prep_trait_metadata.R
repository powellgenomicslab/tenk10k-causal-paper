# prepare trait metadata for workflow
# main input from Anne: https://github.com/powellgenomicslab/TenK10K_SMR/blob/main/metadata/TenK10K_SMR_Metadata.xlsx
# update with additional traits from Albert Henry

# load libraries
library(tidyverse)
library(readxl)
library(fs)
library(janitor)

df_disease <- read_excel("TenK10K_SMR/metadata/TenK10K_SMR_Metadata.xlsx", sheet = "Disease Phenotypes") %>% 
    clean_names() %>% 
    select(trait_id, study, doi, n_case = case, n_control = control) %>%
    mutate(n_total = n_case + n_control,
           n_eff  = floor(4 / (1/n_case + 1/n_control)))

df_bio <- read_excel("TenK10K_SMR/metadata/TenK10K_SMR_Metadata.xlsx", sheet = "Biological measurements") %>% 
    filter(!is.na(trait_id)) %>% 
    clean_names() %>% 
    select(trait_id, study, doi, n_total = sample_size, n_eff = sample_size)

df_otp <- read_tsv("resources/misc/trait_metadata_otp.tsv") %>% 
    mutate(query_id = ifelse(is.na(otp_id), efo_id, otp_id))

# update
df_update <- read_tsv("resources/misc/trait_metadata_n_manual_update.tsv") %>% 
    select(trait_id, study, doi, n_case, n_control, n_total, n_eff)

df_traits <- bind_rows(df_disease, df_bio) %>% 
    filter(!trait_id %in% df_update$trait_id) %>% 
    bind_rows(df_update)

trait_exclude <- c("eo_p", "baso_p", "lymph_p")

df_metadata <- df_otp %>% 
    filter(!trait_id %in% trait_exclude) %>%
    select(trait_id, name, description, supercategory, category, subcategory, efo_id, otp_id, query_id) %>%
    left_join(df_traits)

# write
write_tsv(df_metadata, snakemake@output[[1]])