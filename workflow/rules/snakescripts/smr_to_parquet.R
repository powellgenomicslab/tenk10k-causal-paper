library(arrow)

read_csv_arrow(snakemake@input[[1]]) |>
    write_parquet(snakemake@output[[1]], compression = "gzip")