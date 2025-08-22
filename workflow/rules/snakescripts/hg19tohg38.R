## Purpose: Liftover SNP coordinates from hg19 to hg38 using chain file
## Input: Dataframe with SNP coordinates, chain file
## Output: Dataframe with hg38 coordinates

library(data.table)
library(tidyverse)
library(GenomicRanges)

liftover2hg38 <- function(x, chain_file) {
    #' Liftover SNP coordinates from hg19 to hg38
    #'
    #' @param x Data.table or data.frame with columns: CHROM (chromosome), start (hg19 position), snpid (SNP identifier)
    #' @param chain_file Path to UCSC chain file for hg19-to-hg38 liftover
    #' @return Data.table with columns: chr (chromosome, numeric), pos_hg38 (hg38 position), snpid (SNP identifier)
    #' @details This function takes SNP coordinates in hg19, uses the UCSC chain file to convert them to hg38 coordinates, and returns a table suitable for downstream GWAS or annotation workflows. Input columns must be named as specified, or will be auto-detected if missing. Uses rtracklayer and GenomicRanges for liftover.
    chain <- rtracklayer::import.chain(chain_file)

    # Liftover to hg38 position
    # x should be a dataframe formatted in the following format:
    # data.table(CHROM = chromosome, POS = position (hg19), snpid = SNP identifier to link back to the original gwas)
    x <- as.data.table(x)
    x[, end := start] # Declare end variable
    x[, strand := "*"] # Declare strand variable

    # Reorder columns
    x <- x[, .(seqnames, start, end, strand, snpid)]
    coord_df <- as.data.frame(x)
    coord_granges <- GRanges(coord_df)
    genome(coord_granges) <- "hg19"
    x <- x[, .SD, .SDcols = c("seqnames", "start", "end", "strand", "snpid")]
    # Liftover to hg19
    seqlevelsStyle(coord_granges) <- "UCSC" # necessary
    coord_granges_hg38 <- rtracklayer::liftOver(coord_granges, chain)
    coord_granges_hg38 <- unlist(coord_granges_hg38)
    genome(coord_granges_hg38) <- "hg38"

    # Back to a data table
    hg38_coords <- as.data.table(data.frame(coord_granges_hg38))
    hg38_coords <- hg38_coords[, seqnames := gsub("chr", "", seqnames)]
    hg38_coords <- hg38_coords[, .(seqnames, start, snpid)]
    colnames(hg38_coords) <- c("chr", "pos_hg38", "snpid")
    hg38_coords[, chr := as.numeric(chr)]

    return(hg38_coords)
}

hg38_coords <- hg38_coords[, .SD, .SDcols = c("seqnames", "start", "snpid")]
