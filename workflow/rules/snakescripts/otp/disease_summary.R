## Purpose: Summarize Open Targets target-disease evidence
## Input: Open Targets disease dataset directory
## Output: Disease summary TSV file

library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params


ds_disease <- open_dataset(INPUT$otp_disease_dir)

df_disease <- ds_disease %>%
    collect()

# get summaries of target-diseae pairs
df_areas <- df_disease %>% 
    filter(ontology$isTherapeuticArea == TRUE) %>% 
    select(therapeuticAreas = id,  therapeuticAreaName = name)

df_disease_summary <- df_disease %>%
    filter(ontology$isTherapeuticArea == FALSE) %>%
    select(id, name, therapeuticAreas) %>% 
    unnest(c(therapeuticAreas)) %>% 
    left_join(df_areas)

fs::dir_create(PARAMS$output_dir)
write_tsv(df_disease_summary, OUTPUT[[1]])
