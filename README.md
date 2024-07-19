---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# lazyseurat

<!-- badges: start -->
<!-- badges: end -->

The goal of lazyseurat is to allow import and export Seurat objects to and from
    file-backed storage. 
    Note that LazySeurat is not meant to be used during analysis, 
    but rather for data visualisation on research dissemination
    platforms, where it can save computational costs.

## Installation

You can install the development version of lazyseurat from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jeskowagner/LazySeurat")
```

## Example

Convert a Seurat object to allow lazy loading:


``` r
library(lazyseurat)
library(SeuratData)
#> Error in library(SeuratData): there is no package called 'SeuratData'
# Get example data
InstallData("pbmc3k")
#> Error in InstallData("pbmc3k"): could not find function "InstallData"
data("pbmc3k")
#> Warning in data("pbmc3k"): data set 'pbmc3k' not found

# Export Seurat object to feather for lazy loading
write_seurat_to_feather("lazyseurat")
#> Error in obj[[]]: missing subscript
```

