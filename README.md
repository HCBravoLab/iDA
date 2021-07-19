
# iDA

<!-- badges: start -->
<!-- badges: end -->

This package provides a method for dimensionality reduction
of data to better find discrete latent structure as an improvement over PCA. Developed for the analysis of gene expression where either genetic heterogeneity in the subject population, or cell type in scRNA-seq determines that latent structure. Embeddings are linear transformations so appropriate to use in downstream analyses.

A manuscript describing this method is forthcoming. The package vignette illustrates use of the `iDA` function.

## Installation

The package is available only on github, easiest way to install is using the `remotes` package:

```r
install.packages("remotes")
remotes::install_github("hcbravolab/iDA")
```

