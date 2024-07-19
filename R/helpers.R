#' create_parent_dir
#'
#' @param path Path to create, will only create if it does not exist
#'
#' @export
create_parent_dir <- function(path){

    if(!dir.exists(dirname(path))){
        dir.create(dirname(path), recursive = TRUE)
    }
}

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
get_dense_layer <- function(obj, assay=NULL, layer=NULL){
    counts_raw = LayerData(obj, assay = assay, layer = layer)
    # suppress warnings of large memory allocation, to be expected
    counts_dt = suppressWarnings(transpose(as.data.table(counts_raw)))
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
get_embedding_names <- function(obj){
    embeddings = c()
    for(n in names(obj)){
        if("DimReduc" %in% class(obj[[n]])){
            embeddings = c(embeddings, n)
        }
    }
    embeddings
}


#' subset_data
#' Subset data based on values. Convenience wrapper.
#'
#' @param data The data.frame to subset
#' @param subset_col Which column to use for filtering.
#' @param subset_values Which values to keep.
#'
#' @export
subset_data <- function(data, subset_col, subset_values){
  data[data[[subset_col]] %in% subset_values, ]
}
