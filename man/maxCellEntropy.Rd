% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellEntropy.R
\name{maxCellEntropy}
\alias{maxCellEntropy}
\title{Calculate the max entropy for the cells when given the same counts.}
\usage{
maxCellEntropy(x)
}
\arguments{
\item{x}{A sparse matrix which columns represent cells and rows represent features.}
}
\value{
A vector of the max cell entropy.
}
\description{
Calculate the max entropy for the cells when given the same counts.
}
\examples{
x <- simSimpleCounts()
maxCE <- maxCellEntropy(x)
head(maxCE)
}
