#' table_in_db
#' @description Check if table exists in database.
#'
#' @param db_file Path to database file.
#' @param schema Schema name to check
#' @param table Table name to check
#'
#' @importFrom DBI Id dbListTables
check_table_in_db <- function(db_file, schema, table) {
  con <- withr::local_db_connection(get_connection(db_file))
  if (!table %in% dbListTables(con, Id(schema = schema, table = table))) {
    stop(paste0("Table ", table, " not found in database"))
  }
}


#' Get Schemas in Database
#'
#' @description This function retrieves the names of all schemas
#' in the connected database, excluding the default schemas.
#'
#' @param db_file Path to database file.
#'
#' @return A data frame containing the names of the schemas in the database.
#'
#' @export
get_schemas_in_db <- function(db_file) {
  con <- withr::local_db_connection(get_connection(db_file))
  schemas <- DBI::dbGetQuery(con, "SELECT DISTINCT table_schema FROM information_schema.tables")
  if (nrow(schemas) > 0) {
    return(schemas$table_schema)
  }
  return(NULL)
}

#' Get Tables per Schema in Database
#'
#' @description This function retrieves the names of all tables
#' in the connected database, listing them per schema.
#'
#' @param db_file Path to database file.
#'
#' @return A data frame containing the names of the tables and schemas in the database.
#'
#' @export
get_tables_per_schema_in_db <- function(db_file) {
  con <- withr::local_db_connection(get_connection(db_file))
  DBI::dbGetQuery(con, "SELECT table_schema, table_name FROM information_schema.tables")
}

#' Get Tables for a Specific Schema
#'
#' @description This function retrieves the names of all tables
#' in the connected database, selecting only tables matching a schema.
#'
#' @param db_file Path to database file.
#' @param schema Schema to find tables in
#'
#' @return Character vector with table names, if any.
#'
#' @export
get_tables_in_schema <- function(db_file, schema) {
  get_tables_per_schema_in_db(db_file) %>%
    filter(table_schema == schema) %>%
    pull(table_name) %>%
    unique()
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
#' @param db_file Path to database file.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_metadata <- function(db_file) {
  con <- withr::local_db_connection(get_connection(db_file))
  tbl(con, Id(schema = "metadata", table = "metadata")) %>% return_tibble()
}

#' read_layer
#' @description Read count data of Seurat object.
#'
#' @param db_file Path to database_file
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_layer <- function(db_file, layer = "counts", col_select = NULL) {
  con <- withr::local_db_connection(get_connection(db_file))
  tbl(con, Id(schema = "layer", table = layer)) %>% return_tibble(col_select)
}

#' read_embedding
#' @description Read embedding data of Seurat object.
#'
#' @param db_file Path to database file.
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return tibble of read data
#' @export
#' @importFrom dplyr tbl
#' @importFrom DBI Id
read_embedding <- function(db_file, layer = "umap", col_select = NULL) {
  con <- withr::local_db_connection(get_connection(db_file))
  tbl(con, Id(schema = "embedding", table = layer)) %>%
    return_tibble(col_select)
}

#' read_data_with_meta
#' @description Read count data of Seurat object and merge it with metadata.
#'
#' @param db_file Path to database file.
#' @param what Whether to read count data (`layer`) or an `embedding`.
#' @param name Which layer of data to read (e.g. `counts`).
#' @param col_select Optional vector of columns to read of layer/embedding.
#'
#' @return data.frame of read data
#' @export
#' @importFrom dplyr bind_cols %>%
#' @importFrom duckdb duckdb
read_data_with_meta <- function(db_file,
                                what = "layer",
                                name = "counts",
                                col_select = NULL) {
  check_table_in_db(db_file, schema = "metadata", table = "metadata")
  check_table_in_db(db_file, schema = what, table = name)
  reader <- switch(what,
    layer = read_layer,
    embedding = read_embedding
  )
  read_metadata(db_file) %>% bind_cols(reader(db_file, name, col_select))
}

#' get_connection
#' @description Get connection to database.
#'
#' @param db_file Path to database
#' @param read_only Whether to open connection in read-only mode.
#'
#' @return Connection to database
#' @export
#' @importFrom duckdb duckdb
#' @importFrom DBI dbConnect
get_connection <- function(db_file, read_only = TRUE) {
  dbConnect(duckdb(), db_file, read_only = read_only)
}
