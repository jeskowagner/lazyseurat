library(arrow)

#' read_metadata
#' Read metadata of Seurat object.
#'
#' @param datadir Path to parent directory of feather files.
#' @param col_select Optional vector of columns to read.
#'
#' @return data.frame of read data
#' @export
#' @importFrom arrow read_feather
read_metadata <- function(datadir = ".", col_select = NULL) {
  file <- file.path(datadir, "metadata", "metadata.feather")
  read_feather(file, col_select = col_select)
}

#' read_layer
#' Read count data of Seurat object.
#'
#' @param datadir Path to parent directory of feather files.
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return data.frame of read data
#' @export
#' @importFrom arrow read_feather
read_layer <- function(datadir = ".", layer = "counts", col_select = NULL) {
  file <- file.path(datadir, "layer", paste0(layer, ".feather"))
  read_feather(file, col_select = col_select)
}

#' read_embedding
#' Read embedding data of Seurat object.
#'
#' @param datadir Path to parent directory of feather files.
#' @param layer Which layer to read.
#' @param col_select Optional vector of columns to read.
#'
#' @return data.frame of read data
#' @export
#' @importFrom arrow read_feather
read_embedding <- function(datadir = ".", layer = "umap", col_select = NULL) {
  file <- file.path(datadir, "embedding", paste0(layer, ".feather"))
  read_feather(file, col_select = col_select)
}

#' read_data_with_meta
#' Read count data of Seurat object and merge it with metadata.
#'
#' @param datadir Path to parent directory of feather files.
#' @param what Whether to read count data (`layer`) or an `embedding`.
#' @param name Which layer of data to read (e.g. `counts`).
#' @param col_select Optional vector of columns to read of layer/embedding.
#'
#' @return data.frame of read data
#' @export
#' @importFrom arrow read_feather
read_data_with_meta <- function(datadir = ".",
                                what = "layer",
                                name = "counts",
                                col_select = NULL) {
  reader <- switch(what,
    "layer" = read_layer,
    "embedding" = read_embedding,
    stop("`what` must be 'layer' or 'embedding'")
  )

  metadata <- read_metadata(datadir = datadir)
  data <- reader(datadir = datadir, layer = name, col_select = col_select)
  data <- cbind(metadata, data)
  data
}
