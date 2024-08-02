#' write_seurat_counts
#' @description Write out all data stored in layers, like counts and scaled data.
#'
#' @param obj Seurat object.
#' @param con Connection to database.
#' @param layers Which layers to export. Leave empty for all.
#'
#' @export
#' @importFrom SeuratObject Layers
#' @importFrom DBI Id dbSendQuery dbWriteTable
write_seurat_counts <- function(obj, con, layers = NULL) {
  if (is.null(layers)) {
    layers <- Layers(obj)
  }

  dbSendQuery(con, "CREATE SCHEMA IF NOT EXISTS layer")
  for (l in layers) {
    counts <- get_dense_layer(obj, layer = l)
    dbWriteTable(con,
      name = Id(schema = "layer", table = l),
      value = counts,
      overwrite = TRUE
    )
  }
}

#' write_seurat_embeddings
#' @description Write out all data stored in embeddings, like PCA and UMAP coordinates
#'
#' @param obj Seurat object.
#' @param con Connection to database.
#' @param layers Which layers to export. Leave empty for all.
#'
#' @export
#' @importFrom SeuratObject Embeddings
#' @importFrom DBI Id dbSendQuery dbWriteTable
write_seurat_embeddings <- function(obj, con = ".", layers = NULL) {
  if (is.null(layers)) {
    layers <- get_embedding_names(obj)
  }

  dbSendQuery(con, "CREATE SCHEMA IF NOT EXISTS embedding")
  for (e in layers) {
    coordinates <- as.data.frame(Embeddings(obj, reduction = e))
    dbWriteTable(con,
      name = Id(schema = "embedding", table = e),
      value = coordinates,
      overwrite = TRUE
    )
  }
}

#' write_seurat_metadata
#' @description Write out all data stored in metadata of Seurat object.
#'
#' @param obj Seurat object.
#' @param con Connection to database.
#'
#' @export
#' @importFrom janitor clean_names
#' @importFrom DBI Id dbSendQuery dbWriteTable
write_seurat_metadata <- function(obj, con = ".") {
  metadata <- obj[[]]
  barcode <- rownames(metadata)
  metadata <- cbind(barcode, metadata)
  rownames(metadata) <- NULL
  metadata <- clean_names(metadata)
  dbSendQuery(con, "CREATE SCHEMA IF NOT EXISTS metadata")
  dbWriteTable(con,
    name = Id(schema = "metadata", table = "metadata"),
    value = metadata,
    overwrite = TRUE
  )
}

#' write_seurat_to_db
#' @description Write out all data stored in Seurat object.
#'
#' @param obj Seurat object.
#' @param db Path to database
#'
#' @return Vector of written file paths.
#' @export
#' @importFrom duckdb duckdb
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom withr defer local_db_connection
write_seurat_to_db <- function(obj, db = "seurat.duckdb") {
  con <- local_db_connection(duckdb::dbConnect(duckdb::duckdb(), db, read_only = FALSE))
  defer(suppressWarnings(DBI::dbDisconnect(con, shutdown = TRUE)))
  write_seurat_metadata(obj, con)
  write_seurat_counts(obj, con)
  write_seurat_embeddings(obj, con)
  suppressWarnings(dbDisconnect(con, shutdown = TRUE))
  suppressWarnings(gc())
  invisible(db)
}
