% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{calculateEvaluation}
\alias{calculateEvaluation}
\title{calculateEvaluation(objeto.seurat, mode = c("GTEM", "procrustes", "RVSSimageJ"),
listaCoordenadasNEW, listaCoordenadas,
patientType = c('unique','multiple'))}
\usage{
calculateEvaluation(
  objeto.seurat,
  mode = c("GTEM", "procrustes", "RVSSimageJ"),
  listaCoordenadasNEW = NULL,
  listaCoordenadas = NULL,
  patientType = c("unique", "multiple")
)
}
\arguments{
\item{objeto.seurat}{Seurat object containing images and coordinates.}

\item{mode}{Evaluation mode: one of "GTEM", "procrustes", or "RVSSimageJ".}

\item{listaCoordenadasNEW}{List of new coordinates (optional, used mainly if mode != "RVSSimageJ").}

\item{listaCoordenadas}{List of original coordinates (optional, used mainly if mode != "RVSSimageJ").}

\item{patientType}{Patient type, affecting region size ("unique" or "multiple").}
}
\value{
Evaluation results as a list.
}
\description{
Evaluate the alignment of images using various metrics (MSE, SSIM, etc.) for both raw and transformed images.
}
