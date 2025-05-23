% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{EuclDist}
\alias{EuclDist}
\title{EuclDist(listaCoordenadas, nIm)}
\usage{
EuclDist(listaCoordenadas, nIm)
}
\arguments{
\item{listaCoordenadas}{A list where each element contains a list with elements $x and $y representing coordinates.}

\item{nIm}{A list or vector of length 2 with indices of the images to compare.}
}
\value{
Numeric value representing the Euclidean distance between the two coordinate sets.
}
\description{
Calculate Euclidean distance between two sets of coordinates.
}
