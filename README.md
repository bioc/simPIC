# simPIC

## Overview
<img src='vignettes/logo/Logo.png' align="right" height="139" />

simPIC is an R package for simple simulation of single-cell Assay for 
Transposase Accessible Chromatin sequencing (scATAC-seq) data. 
simPIC provides an easy to use interface for:

* estimating simulation parameters
* Objects for storing those parameters
* simulating counts using those parameters

## News

**Version 1.5.3 (Development Version)**

- Major updates to `simPICsimulate` function to allow simulating multiple 
cell-types and batch effects.

For full change logs, please refer to the [NEWS file](https://github.com/sagrikachugh/simPIC/blob/devel/NEWS.md).

## Installation

The package can be installed from Bioconductor using the following code

```r
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("simPIC")
```
For development version

```r
BiocManager::install(
    "sagrikachugh/simPIC",
    dependencies = TRUE,
    build_vignettes = TRUE
)
```

## Getting started

To get started, check out the vignette for a quick start and detailed look into
simPIC.

```r
library(simPIC)
browseVignettes("simPIC")
```


