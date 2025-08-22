library(fs)
library(tidyverse)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

calc_q <- function(p, ...) {
    possibly(qvalue, qvalue_truncp(p))(p, ...) %>% .$qvalues
}


df_cell_type_stats <- map_df(
    INPUT[['cell_type_stats']] %>% 
        set_names(., str_remove(basename(.), "\\.cell_type_stats\\.tsv$")),
    ~read_tsv(.x) %>% 
        rename(cell_type = group) %>% 
        mutate(qval_assoc_mcp = calc_q(assoc_mcp),
               bh_assoc_mcp = p.adjust(assoc_mcp, method = "BH")),
    .id = "phenotype"
)

df_cell_type_top <-  map_df(
    INPUT[['cell_type_top']] %>% 
        set_names(., str_remove(basename(.), "\\.cell_type_top\\.tsv$")),
    ~read_tsv(.x),
    .id = "phenotype"
)

write_tsv(df_cell_type_stats, OUTPUT[['cell_type_stats']])
write_tsv(df_cell_type_top, OUTPUT[['cell_type_top']])
