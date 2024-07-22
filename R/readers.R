library(duckdb)
library(DBI)

#' table_in_db
#' Check if table exists in database.
#'
#' @param con Connection to database
#' @param table Table name to check
check_table_in_db <- function(con, schema, table) {
  if (!table %in% dbListTables(con, Id(schema = schema, table = table))) {
    stop(paste0("Table ", table, " not found in database"))
  }
}

#' return_tibble
#' Collect and return data as tibble.
#'
#' @param table Table to collect
#' @param col_select Optional vector of columns to select.
#'
#' @return tibble of collected data
#' @importFrom dplyr collect as_tibble all_of select %>%
return_tibble <- function(table, col_select = NULL) {
  if (!is.null(col_select)) {
    table <- table %>% select(all_of(col_select))
  }
  table %>%
    collect() %>%
    as_tibble()
}

#' read_metadata
#' Read metadata of Seurat object.
#'
#' @param con Connection to database
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
read_metadata <- function(con) {
  tbl(con, Id(schema = "metadata", table = "metadata")) %>% return_tibble()
}

#' read_layer
#' Read count data of Seurat object.
#'
#' @param con Connection to database
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
read_layer <- function(con, layer = "counts", col_select = NULL) {
  tbl(con, Id(schema = "layer", table = layer)) %>% return_tibble(col_select)
}

#' read_embedding
#' Read embedding data of Seurat object.
#'
#' @param con Connection to database
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
read_embedding <- function(con, layer = "umap", col_select = NULL) {
  tbl(con, Id(schema = "embedding", table = layer)) %>%
    return_tibble(col_select)
}

#' read_data_with_meta
#' Read count data of Seurat object and merge it with metadata.
#'
#' @param db Path to database
#' @param what Whether to read count data (`layer`) or an `embedding`.
#' @param name Which layer of data to read (e.g. `counts`).
#' @param col_select Optional vector of columns to read of layer/embedding.
#'
#' @return data.frame of read data
#' @export
#' @importFrom dplyr bind_cols %>%
read_data_with_meta <- function(db,
                                what = "layer",
                                name = "counts",
                                col_select = NULL) {
  con <- dbConnect(duckdb(), db, read_only = TRUE)
  check_table_in_db(con, schema = "metadata", table = "metadata")
  check_table_in_db(con, schema = what, table = name)
  reader <- switch(what,
    layer = read_layer,
    embedding = read_embedding
  )
  res <- read_metadata(con) %>% bind_cols(reader(con, name, col_select))
  dbDisconnect(con)
  return(res)
}
