# RcwlPipelines

The RcwlPipelines package manages a collection of Bioinformatics tools and pipeline recipes based on Rcwl.

## Installation

The package can be installed from Bioconductor (>= 3.9):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RcwlPipelines")
	
# Or from github
BiocManager::install("hubentu/RcwlPipelines")
```

## User Guide

``` r
vignette(package = "RcwlPipelines")
```
