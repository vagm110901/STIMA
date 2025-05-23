% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/publishedMethodsComparison.R
\name{createDeconvolutionListsPython}
\alias{createDeconvolutionListsPython}
\title{deconvolution
createDeconvolutionListsPython(object.seurat, split_by = "name", mode = c("STalign", "PASTE2"))}
\usage{
createDeconvolutionListsPython(
  object.seurat,
  object,
  split_by = "name",
  mode = c("STalign", "PASTE2")
)
}
\arguments{
\item{object.seurat}{A Seurat object after the alignment.
If mode == STalign, it should be the Seurat object after the STalign alignment. {STaligntoSeurat()}
If mode == PASTE2, it should be the Seurat object after the PASTE2 alignment. {PASTE2toSeurat()}}

\item{object}{A Seurat object with the original data.}

\item{split_by}{A character string indicating the metadata field to split the Seurat object by. Default is "name".}

\item{mode}{A character string indicating the deconvolution method to be used. Options are "GTEM", "procrustes", or "RVSSimageJ". Default is "GTEM".}
}
\value{
A list of Seurat objects, each containing a reference and two problematic objects for deconvolution analysis.
}
\description{
Creates a list of Seurat objects for deconvolution analysis.
It uses an objetc created by python (STalign or PASTE2).
}
