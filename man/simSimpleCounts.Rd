% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simSimpleCounts}
\alias{simSimpleCounts}
\title{Simulate simple counts matrix for the empty and cell-containing droplets.}
\usage{
simSimpleCounts(
  total_gene = 30000,
  nempty = 20000,
  nlarge = 2000,
  nsmall = 200,
  empty_prof = NULL,
  empty_ngene_rate = 0.05,
  empty_rate = 1/100,
  large_prof = NULL,
  large_ngene_rate = 0.8,
  large_shape = 6,
  large_scale = 500,
  small_prof = NULL,
  small_ngene_rate = 0.6,
  small_shape = 10,
  small_scale = 100,
  remove_zero_drop = TRUE,
  remove_zero_feature = TRUE,
  seed = 0
)
}
\arguments{
\item{total_gene}{Total gene number for all the simulated counts.}

\item{nempty}{Empty, large cell and small cell droplets number. If \code{remove_zero_drop} is \code{TRUE}, the final number may be smaller beacuse droplets that have all zero-valued counts will be removed.}

\item{nlarge}{Empty, large cell and small cell droplets number. If \code{remove_zero_drop} is \code{TRUE}, the final number may be smaller beacuse droplets that have all zero-valued counts will be removed.}

\item{nsmall}{Empty, large cell and small cell droplets number. If \code{remove_zero_drop} is \code{TRUE}, the final number may be smaller beacuse droplets that have all zero-valued counts will be removed.}

\item{empty_prof}{The overall gene expression profile distribution. If provided, must be the same length with \code{total_gene}.}

\item{empty_ngene_rate}{Rate of total genes expressed in each type of droplets.}

\item{empty_rate}{Rate parameter of exponential distribution for 'Empty'.}

\item{large_prof}{The overall gene expression profile distribution. If provided, must be the same length with \code{total_gene}.}

\item{large_ngene_rate}{Rate of total genes expressed in each type of droplets.}

\item{large_shape}{shape parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale} control the expected mean value of counts in large or small cells.}

\item{large_scale}{scale parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale^2} control the expected variance of counts in large or small cells.}

\item{small_prof}{The overall gene expression profile distribution. If provided, must be the same length with \code{total_gene}.}

\item{small_ngene_rate}{Rate of total genes expressed in each type of droplets.}

\item{small_shape}{shape parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale} control the expected mean value of counts in large or small cells.}

\item{small_scale}{scale parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale^2} control the expected variance of counts in large or small cells.}

\item{remove_zero_drop}{Whether to remove all zero-valued droplets.}

\item{remove_zero_feature}{Whether to remove all zero-valued features.}

\item{seed}{Random seed used in simulation. Default is 0.}
}
\value{
A sparse Matrix of class "dgCMatrix".
}
\description{
Simulate simple counts matrix for the empty and cell-containing droplets.
}
\examples{
counts <- simSimpleCounts()
counts
}
