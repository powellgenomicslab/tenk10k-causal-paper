# misc helper function

library(dplyr)

write_table <- function(df, source_tbl_name,
                        tbl_order = NULL,
                        target_tbl_name = source_tbl_name,
                        target_dir = "tables",
                        cols = NULL,
                        label_xlsx = "metadata/table_column_names.xlsx") {
    if (!is.null(cols)) {
        df <- df %>% select(all_of(cols))
    }
    if (!is.null(tbl_order)) {
        target_tbl_name <- paste0(tbl_order, "-", target_tbl_name)
    }
    labels <- readxl::read_excel(label_xlsx, sheet = source_tbl_name) %>%
        select(label, name) %>%
        deframe()
    df %>%
        select(!!!labels) %>%
        data.table::fwrite(file.path(target_dir, paste0(target_tbl_name, ".tsv")), sep = "\t", na = "NA", quote = FALSE)
}
