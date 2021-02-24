# RcwlPipelines

The RcwlPipelines package manages a collection of Bioinformatics tools and pipeline recipes based on Rcwl. The pre-built and pre-tested tools and pipelines are highly modularized with easy customization to meet different bioinformatics data analysis needs.

## Installation

The package can be installed from Bioconductor (>= 3.9):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RcwlPipelines")
	
# Or from github
BiocManager::install("rworkflow/RcwlPipelines")
```

## Get started

``` r
cwlUpdate()
cwlSearch("STAR")
STAR <- cwlLoad("tl_STAR")
```

## User Guide

``` r
vignette(package = "RcwlPipelines")
```
