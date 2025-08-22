## Purpose: Summarize known drug information from Open Targets
## Input: Known drug dataset directory
## Output: Drug summary TSV file

library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params


ds_drug <- open_dataset(INPUT$otp_drug_dir)

df_drug <- ds_drug %>%
    collect()

# get summaries of target-diseae pairs
df_drug_summary <- df_drug %>%
    mutate(n_trial = map_dbl(urls, length)) %>%
    select(where(~ !is.list(.x)))

fs::dir_create(PARAMS$output_dir)
write_tsv(df_drug_summary, OUTPUT$drug_summary)
