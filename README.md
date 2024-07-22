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

## Example

### Step 1: Obtain some example data

``` r
# Option 1: load in your own dataset
library(Seurat)
obj = readRDS("path/to/your/seurat_object.rds")

# Option 2: use the pbmc3k dataset
# Install SeuratData to get example data
if(!requireNamespace("SeuratData", quietly = TRUE)){
    devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)
# Download example data
if("pbmc3k.SeuratData" %in% rownames(InstalledData())){
    data("pbmc3k")
} else {
    InstallData("pbmc3k")
}
obj = pbmc3k
```

### Step 2: Convert Seurat object to database


``` r
library(lazyseurat)

# Export Seurat object to DuckDB
# Note: this may take a minute
write_seurat_to_db(obj, "seurat.duckdb")
```

### Step 3: Load Seurat data from database


``` r
# Open connection to database
con <- get_connection("seurat.duckdb")

# Read MYC gene expression
read_data_with_meta(con=con, col_select="MYC")
```
