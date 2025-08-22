# quantify open targets target-disease evidence

# Using arrow / dplyr
library(tidyverse)
library(arrow)
library(readxl)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
df_trait_meta <- read_excel(INPUT$trait_metadata)

# helper function to query a dataset

filter_unique <- function(dataset, col, x) {
    x <- discard(unique(x), ~ is.na(.x))
    dataset %>%
        filter({{ col }} %in% x)
}

df_drug <- INPUT$otp_drug_dir %>%
    open_dataset() %>%
    filter_unique(diseaseId, df_trait_meta$query_id) %>%
    collect()

# get summaries of target-diseae pairs
df_drug_summary <- df_drug %>%
    unnest(urls) %>%
    group_by(targetId, diseaseId) %>%
    summarise(max_phase = max(phase, na.rm = TRUE),
              n_drug = length(unique(drugId)),
              mean_phase = mean(phase, na.rm = TRUE),
              n_trial = length(url)) %>%
    ungroup()

fs::dir_create(PARAMS$output_dir)
write_parquet(df_drug, OUTPUT$drug, compression = "gzip")
write_parquet(df_drug_summary, OUTPUT$drug_summary, compression = "gzip")
