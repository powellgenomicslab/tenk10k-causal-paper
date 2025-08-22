## Purpose: Generate gene annotation table from GTF file
## Input: GTF file
## Output: Gene annotation TSV file
library(data.table)
library(dplyr)
library(stringr)

rtracklayer::import(snakemake@input[[1]]) %>% 
  as.data.table() %>% 
  filter(type == "gene") %>% 
  mutate(ensembl_gene_id = stringr::str_remove_all(gene_id, "\\..*"),
         chr = stringr::str_remove(seqnames, "^chr") %>% as.numeric) %>% 
  filter(!is.na(chr)) %>% 
  select(chr, start, end, ensembl_gene_id, hgnc_symbol = gene_name, gene_type) %>% 
  fwrite(snakemake@output[[1]], sep = "\t")