% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellCalling.R
\name{CallCellRangerV3}
\alias{CallCellRangerV3}
\title{CellRanger_v3 method used to recognize cell-containing droplets.}
\usage{
CallCellRangerV3(
  counts,
  recovered_cells = 3000,
  recovered_cells_quantile = 0.99,
  n_candidate_barcodes = 20000,
  n_partitions = 90000,
  min_umis_nonambient = 500,
  min_umi_frac_of_median = 0.01,
  max_adj_pvalue = 0.01
)
}
\arguments{
\item{counts}{A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.}

\item{recovered_cells}{Expected number of recovered cells. Default is 3000.}

\item{recovered_cells_quantile}{Quantile of the top \code{recovered_cells} barcodes by total UMI counts. Default is 0.99.}

\item{n_candidate_barcodes}{Number of additional barcodes to consider after the initial cell calling. Default is 20000.}

\item{n_partitions}{Number of partitions (max number of barcodes to consider for ambient estimation). Default is 90000.}

\item{min_umis_nonambient}{Minimum number of UMIS per barcode to consider after the initial cell calling. Default is 500.}

\item{min_umi_frac_of_median}{Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling. Default is 0.01.}

\item{max_adj_pvalue}{Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling. Default is 0.01.}
}
\value{
A \code{DataFrame} containing the classification column named 'CellRangerV3Class'.
}
\description{
CellRanger_v3 method used to recognize cell-containing droplets.
}
\examples{
counts <- simSimpleCounts()
result <- CallCellRangerV3(counts)
head(result)
}
