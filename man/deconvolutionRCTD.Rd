% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolution.R
\name{deconvolutionRCTD}
\alias{deconvolutionRCTD}
\title{deconvolutionRCTD()}
\usage{
deconvolutionRCTD(
  object.list,
  reference,
  im,
  mode = c("GTEM", "procrustes", "RVSSimageJ")
)
}
\arguments{
\item{object.list}{A list of Seurat objects containing the data to be deconvoluted.}

\item{reference}{Object with the single-cell / single-nucelus reference for the deconvolution.}

\item{im}{A numeric value indicating the image number.}

\item{mode}{A character string indicating the deconvolution method to be used. Options are "GTEM", "procrustes", or "RVSSimageJ". Default is "GTEM".}
}
\value{
A list of Seurat objects with deconvolution results added as assays.
}
\description{
Performs deconvolution using the RCTD package for a list of reference-NOTalign-YESalign.
}
