
# STIMA

<!-- badges: start -->

<!-- badges: end -->

STIMA (Spatial Transcriptomics Image-based Methods for Alignment), is
designed to align two or more ST slices or samples, enabling the
comparison and analysis of gene expression within the same regions.
STIMA performs alignment in a pairwise comparison manner, considering
one slice as a reference, which includes both the tissue microscope
image and the spatial spot matrix of gene expression. However, STIMA
relies exclusively on image data for alignment without incorporating
gene expression data, thereby preserving the independence of
transcriptomic information across samples. STIMA includes three distinct
alignment approaches—Geometric Transformation Estimation Model (GTEM),
Procrustes Transformation and the ImageJ plugin Register Virtual Stack
Slices (RVSS-ImageJ).

## Installation

You can install the development version of STIMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("vagm110901/STIMA")
```
