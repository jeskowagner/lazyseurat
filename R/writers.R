#' write_seurat_counts
#' Write out all data stored in layers, like counts and scaled data.
#'
#' @param obj Seurat object.
#' @param outdir Path to output directory.
#' @param layers Which layers to export. Leave empty for all.
#'
#' @return Vector of written file paths.
#' @export
#' @importFrom arrow write_feather
#' @importFrom SeuratObject Layers
write_seurat_counts <- function(obj, outdir=".", layers=NULL){

    if(is.null(layers)){
        layers = Layers(obj)
    }

    out_paths = c()
    for(l in layers){
        out_paths = c(out_paths, file.path(outdir, "layer", paste0(l, ".feather")))
        create_parent_dir(out_paths[length(out_paths)])
        counts = get_dense_layer(obj, layer=l)
        write_feather(counts, out_paths[length(out_paths)])
        rm(counts)
        gc()
    }
    gc()
    invisible(out_paths)
}

#' write_seurat_embeddings
#' Write out all data stored in embeddings, like PCA and UMAP coordinates
#'
#' @param obj Seurat object.
#' @param outdir Path to output directory.
#' @param layers Which layers to export. Leave empty for all.
#'
#' @return Vector of written file paths.
#' @export
#' @importFrom data.table as.data.table
#' @importFrom arrow write_feather
#' @importFrom SeuratObject Embeddings
write_seurat_embeddings <- function(obj, outdir=".", layers=NULL){
    if(is.null(layers)){
      layers = get_embedding_names(obj)
    }

    out_paths = c()
    for(e in layers){
        out_paths = c(out_paths, file.path(outdir, "embedding", paste0(e, ".feather")))
        create_parent_dir(out_paths[length(out_paths)])
        coordinates = as.data.table(Embeddings(obj, reduction = e))
        write_feather(coordinates, out_paths[length(out_paths)])
        rm(coordinates)
        gc()
    }
    gc()
    invisible(out_paths)
}

#' write_seurat_metadata
#' Write out all data stored in metadata of Seurat object.
#'
#' @param obj Seurat object.
#' @param outdir Path to output directory.
#'
#' @return Vector of written file paths.
#' @export
#' @importFrom arrow write_feather
write_seurat_metadata <- function(obj, outdir="."){
    out_path = file.path(outdir, "metadata", "metadata.feather")
    create_parent_dir(out_path)
    metadata = obj[[]]
    barcode = rownames(metadata)
    metadata = cbind(barcode, metadata)
    write_feather(metadata, out_path)
    invisible(out_path)
}

#' write_seurat_to_feather
#' Write out all data stored in Seurat object.
#'
#' @param obj Seurat object.
#' @param outdir Path to output directory.
#'
#' @return Vector of written file paths.
#' @export
write_seurat_to_feather <- function(obj, outdir="."){
    meta_files = write_seurat_metadata(obj, outdir)
    counts_files = write_seurat_counts(obj, outdir)
    embedding_files = write_seurat_embeddings(obj, outdir)
    invisible(return(c(meta_files, counts_files, embedding_files)))
}
