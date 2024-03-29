% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simComplexCounts}
\alias{simComplexCounts}
\title{Simulate a complex counts matrix including different types of empty and cell-containing droplets.}
\usage{
simComplexCounts(
  total_gene = 30000,
  disturbance = 0.2,
  nempty = 50000,
  nlarge = 5000,
  nsmall = 500,
  empty_type = 2,
  large_type = 10,
  small_type = 2,
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
  large_frag = TRUE,
  large_frag_gene = 1:50,
  large_frag_prop = 0.5,
  small_frag = TRUE,
  small_frag_gene = 1:50,
  small_frag_prop = 0.5,
  remove_zero_drop = TRUE,
  remove_zero_feature = TRUE,
  seed = 0
)
}
\arguments{
\item{total_gene}{Total gene number for all the simulated counts.}

\item{disturbance}{A numeric value used as a weight of standard deviations when sample different distribution parameters for each cell/empty type from defined global parameters. Default is 0.2.}

\item{nempty, nlarge, nsmall}{Empty, large cell and small cell droplets number. If \code{remove_zero_drop} is \code{TRUE}, the final number may be smaller beacuse droplets that have all zero-valued counts will be removed.}

\item{empty_type, large_type, small_type}{Total number of types for \code{nempty,nlarge,nsmall}.}

\item{empty_prof, large_prof, small_prof}{The overall gene expression profile distribution. If provided, must be the same length with \code{total_gene}.}

\item{empty_ngene_rate, large_ngene_rate, small_ngene_rate}{Rate of total genes expressed in each type of droplets.}

\item{empty_rate}{Rate parameter of exponential distribution for 'Empty'.}

\item{large_shape, small_shape}{shape parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale} control the expected mean value of counts in large or small cells.}

\item{large_scale, small_scale}{scale parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale^2} control the expected variance of counts in large or small cells.}

\item{large_frag, small_frag}{Whether simulate cell fragments from the large or small cells. Default is TRUE.}

\item{large_frag_gene, small_frag_gene}{Indices of cell fragment gene in profile. Default is 1:100.}

\item{large_frag_prop, small_frag_prop}{Proportion of the cell fragment gene counts. Default is 0.5.}

\item{remove_zero_drop}{Whether to remove all zero-valued droplets.}

\item{remove_zero_feature}{Whether to remove all zero-valued features.}

\item{seed}{Random seed used in simulation. Default is 0.}
}
\value{
A sparse Matrix of class "dgCMatrix".
}
\description{
Simulation of a complex single-cell sequencing dataset.
}
\examples{
counts <- simComplexCounts()
counts
}
