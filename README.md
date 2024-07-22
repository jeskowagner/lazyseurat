<!-- README.md is generated from README.Rmd. Please edit that file -->



# lazyseurat

<!-- badges: start -->
[![R-CMD-check](https://github.com/jeskowagner/LazySeurat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jeskowagner/LazySeurat/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of lazyseurat is to allow import and export Seurat objects to and from
    file-backed storage.
    Note that lazyseurat is not meant to be used during analysis,
    but rather for data visualisation on research dissemination
    platforms, where it can save computational costs.

## Installation

You can install the development version of lazyseurat from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jeskowagner/lazyseurat")
```

<!-- ## Example -->

<!-- Convert a Seurat object to allow lazy loading: -->

<!-- ```{r export_seurat} -->
<!-- library(lazyseurat) -->

<!-- # Install SeuratData to get example data -->
<!-- if(!requireNamespace("SeuratData", quietly = TRUE)){ -->
<!--     devtools::install_github('satijalab/seurat-data') -->
<!-- } -->

<!-- library(SeuratData) -->

<!-- # Download example data -->
<!-- if("pbmc3k.SeuratData" %in% rownames(InstalledData())){ -->
<!--     data("pbmc3k") -->
<!-- } else { -->
<!--     InstallData("pbmc3k") -->
<!-- } -->

<!-- # Export Seurat object to feather for lazy loading -->
<!-- write_seurat_to_feather(pmbc3k, "lazyseurat") -->
<!-- ``` -->
