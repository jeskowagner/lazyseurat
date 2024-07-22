#' get_dense_layer
#' Convert a layer in Seurat object to dense data.table, with rows being cells.
#'
#' @param obj Seurat object
#' @param assay Which assay to use (e.g. "RNA")
#' @param layer Which layer to use (e.g. "counts")
#'
#' @return `data.table` with dense layer information.
#'
#' @export
#' @importFrom data.table transpose setnames as.data.table
#' @importFrom SeuratObject LayerData
get_dense_layer <- function(obj, assay = NULL, layer = NULL) {
  require(SeuratObject)
  counts_raw <- LayerData(obj, assay = assay, layer = layer)
  # suppress warnings of large memory allocation, to be expected
  counts_dt <- suppressWarnings(transpose(as.data.table(counts_raw)))
  setnames(counts_dt, rownames(counts_raw))
  counts_dt
}

#' get_embedding_names
#' Get names of dimensionality reduction embeddings in Seurat object.
#'
#' @param obj Seurat object
#'
#' @return Vector of the embedding names.
#'
#' @export
get_embedding_names <- function(obj) {
  embeddings <- c()
  for (n in names(obj)) {
    if ("DimReduc" %in% class(obj[[n]])) {
      embeddings <- c(embeddings, n)
    }
  }
  embeddings
}

#' check_if_duplicated_columns
#' Check if there are duplicated columns in a data.frame.
#'
#' @param df The data.frame to check
#'
#' @return Vector with names of duplicated columns.
#' @export
#' @importFrom stringr str_to_lower
which_duplicated_columns <- function(df) {
  cols <- str_to_lower(colnames(df))
  dups <- duplicated(cols) | duplicated(cols, fromLast = TRUE)
  ifelse(any(dups), cols[dups], NULL)
}
