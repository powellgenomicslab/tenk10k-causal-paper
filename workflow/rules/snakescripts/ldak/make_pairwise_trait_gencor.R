library(tidyverse)
library(fs)

INPUT <- snakemake@input[[1]]
OUTPUT <- snakemake@output[[1]]
PARAMS <- snakemake@params
traits <- read_lines(INPUT)

# Create a data frame with all combinations of traits
trait_pairs <- expand_grid(trait1 = traits, trait2 = traits) %>%
    filter(trait1 != trait2) %>%
    mutate(pairs = paste(pmin(trait1, trait2), pmax(trait1, trait2), sep = ".")) %>%
    distinct(pairs, .keep_all = TRUE) %>%
    mutate(ldak_cor = path(PARAMS$ldak_prefix,
                         paste(trait1, trait2, "cors", sep = ".")))

dir_create(dirname(OUTPUT))
write_lines(trait_pairs$ldak_cor, OUTPUT)