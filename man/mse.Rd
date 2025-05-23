% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{mse}
\alias{mse}
\title{mse(image1, image2)}
\usage{
mse(image1, image2)
}
\arguments{
\item{image1}{A 3D array representing the first image (height x width x channels).}

\item{image2}{A 3D array representing the second image (same dimensions as image1).}
}
\value{
Numeric value representing the average MSE across RGB channels.
}
\description{
Calculate Mean Squared Error (MSE) between two images by comparing pixel values.
}
