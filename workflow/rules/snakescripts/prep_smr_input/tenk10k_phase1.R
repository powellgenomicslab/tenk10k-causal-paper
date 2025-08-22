## Purpose: Calculate p-threshold and generate probe list for SMR per celltype

library(data.table)
library(fs)
library(tidyverse)

INPUT <- "brenner/tenk10k_phase1/common_eqtl.tsv"
OUTPUT <- snakemake@output

df <- fread(INPUT)

cells <- unique(df$celltype)

process <- function(c) {
    dir_create(path(OUTPUT, c))
    thresh <- df[df$celltype == c, ][, max(Pvalue)]
    probe_list <- df[df$celltype == c, ][, unique(Gene)]
    ## Output is written to pthresh.txt per celltype
    readr::write_lines(as.character(thresh), path(OUTPUT, c, "pthresh.txt"))
    readr::write_lines(as.character(probe_list), path(OUTPUT, c, "probe.txt"))
}

walk(cells, process)