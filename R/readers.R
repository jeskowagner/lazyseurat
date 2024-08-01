#' table_in_db
#' @description Check if table exists in database.
#'
#' @param con Connection to database
#' @param schema Schema name to check
#' @param table Table name to check
#'
#' @importFrom DBI Id dbListTables
check_table_in_db <- function(con, schema, table) {
  if (!table %in% dbListTables(con, Id(schema = schema, table = table))) {
    stop(paste0("Table ", table, " not found in database"))
  }
}

#' return_tibble
#' @description Collect and return data as tibble.
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
#' @description Read metadata of Seurat object.
#'
#' @param con Connection to database
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_metadata <- function(con) {
  tbl(con, Id(schema = "metadata", table = "metadata")) %>% return_tibble()
}

#' read_layer
#' @description Read count data of Seurat object.
#'
#' @param con Connection to database
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_layer <- function(con, layer = "counts", col_select = NULL) {
  tbl(con, Id(schema = "layer", table = layer)) %>% return_tibble(col_select)
}

#' read_embedding
#' @description Read embedding data of Seurat object.
#'
#' @param con Connection to database
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_embedding <- function(con, layer = "umap", col_select = NULL) {
  tbl(con, Id(schema = "embedding", table = layer)) %>%
    return_tibble(col_select)
}

#' read_data_with_meta
#' @description Read count data of Seurat object and merge it with metadata.
#'
#' @param con Connection to database
#' @param what Whether to read count data (`layer`) or an `embedding`.
#' @param name Which layer of data to read (e.g. `counts`).
#' @param col_select Optional vector of columns to read of layer/embedding.
#'
#' @return data.frame of read data
#' @export
#' @importFrom dplyr bind_cols %>%
#' @importFrom duckdb duckdb
read_data_with_meta <- function(con,
                                what = "layer",
                                name = "counts",
                                col_select = NULL) {
  check_table_in_db(con, schema = "metadata", table = "metadata")
  check_table_in_db(con, schema = what, table = name)
  reader <- switch(what,
    layer = read_layer,
    embedding = read_embedding
  )
  read_metadata(con) %>% bind_cols(reader(con, name, col_select))
}

#' get_connection
#' @description Get connection to database.
#'
#' @param db Path to database
#' @param read_only Whether to open connection in read-only mode.
#'
#' @return Connection to database
#' @export
#' @importFrom duckdb duckdb
#' @importFrom DBI dbConnect
get_connection <- function(db, read_only = TRUE) {
  dbConnect(duckdb(), db, read_only = read_only)
}
