% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment.R
\name{resultProcrustes}
\alias{resultProcrustes}
\title{resultProcrustes(proc, mirrorx, mirrory, scale)}
\usage{
resultProcrustes(proc, mirrorx, mirrory, scale)
}
\arguments{
\item{proc}{A list containing the results of a Procrustes analysis. Expected to have elements:
\itemize{
\item \code{R}: a rotation matrix,
\item \code{t}: a translation vector,
\item \code{d}: a scaling factor (distance).
}}

\item{mirrorx}{0 or anothe value; whether to apply reflection along the x-axis.}

\item{mirrory}{0 or anothe value; whether to apply reflection along the y-axis.}

\item{scale}{Logical; whether to apply scaling.}
}
\value{
A named numeric vector with the following elements:
\itemize{
\item \code{coseno}: cosine of the rotation angle,
\item \code{seno}: sine of the rotation angle,
\item \code{dx}: translation along the x-axis,
\item \code{dy}: translation along the y-axis,
\item \code{mirrorx}: reflection flag on the x-axis,
\item \code{mirrory}: reflection flag on the y-axis,
\item \code{e}: scaling factor (1 if \code{scale} is FALSE).
}
}
\description{
Extracts the key transformation parameters from a Procrustes analysis result,
including rotation, translation, reflection (mirroring), and scaling.
}
