#' write_seurat_counts
#' @description Write out all data stored in layers, like counts and scaled data.
#'
#' @param obj Seurat object.
#' @param con Connection to database.
#' @param layers Which layers to export. Leave empty for all.
#'
#' @export
#' @importFrom SeuratObject Layers LayerData
#' @importFrom DBI Id dbSendQuery dbWriteTable
#' @importFrom data.table setnames
#' @importFrom dplyr %>% as_tibble
#' @importFrom tibble rownames_to_column
write_seurat_counts <- function(obj, con, layers = NULL) {
  if (is.null(layers)) {
    layers <- Layers(obj)
  }
  cat("Writing counts to database...\n")
  dbSendQuery(con, "CREATE SCHEMA IF NOT EXISTS layer")
  dbSendQuery(con, "CREATE SCHEMA IF NOT EXISTS averages")

  for (l in layers) {
    counts <- get_dense_layer(obj, layer = l)
    averages <- LayerData(obj, layer = l) %>%
      Matrix::rowMeans(.) %>%
      as.data.frame %>%
      setnames("average") %>%
      rownames_to_column("gene")


    dbWriteTable(con,
      name = Id(schema = "layer", table = l),
      value = counts,
      overwrite = TRUE
    )

    dbWriteTable(con,
      name = Id(schema = "averages", table = l),
      value = averages,
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
  cat("Writing embeddings to database...\n")
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
#' @description Write out all data stored in metadata and idents
#' of Seurat object.
#'
#' @param obj Seurat object.
#' @param con Connection to database.
#'
#' @export
#' @importFrom janitor clean_names
#' @importFrom DBI Id dbSendQuery dbWriteTable
#' @importFrom tibble rownames_to_column
write_seurat_metadata <- function(obj, con = ".") {
  metadata <- obj[[]]
  barcode <- rownames(metadata)
  metadata <- cbind(barcode, metadata)
  rownames(metadata) <- NULL
  idents <- rownames_to_column(as.data.frame(Idents(obj)), "barcode")
  colnames(idents) <- c("barcode", "ident")
  metadata <- merge(metadata, idents, by = "barcode")
  metadata <- clean_names(metadata)

  cat("Writing metadata to database...\n")

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
  defer(suppressWarnings(DBI::dbDisconnect(con, shutdown = TRUE)), priority = "last")
  write_seurat_metadata(obj, con)
  write_seurat_counts(obj, con)
  write_seurat_embeddings(obj, con)
  suppressWarnings(dbDisconnect(con, shutdown = TRUE))
  suppressWarnings(gc())
  invisible(db)
}
