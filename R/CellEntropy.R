#' Calculate entropy for each cell according to the feature counts.
#'
#' @param x A sparse matrix which columns represent cells and rows represent features.
#'
#' @return A vector of the cell entropy.
#'
#' @examples
#' x <- simCounts()
#' CE <- CellEntropy(x)
#' head(CE)
#' @importFrom Matrix t colSums
#' @export
CellEntropy <- function(x) {
  if (!class(x) %in% c("dgCMatrix", "dgTMatrix")) {
    stop("'x' must be sparse Matrix of class dgCMatrix or dgTMatrix")
  }
  logp_matrix <- p_matrix <- Matrix::t(Matrix::t(x) / Matrix::colSums(x))
  logp_matrix@x <- log2(logp_matrix@x)
  e <- p_matrix * logp_matrix
  entropy <- -Matrix::colSums(e, na.rm = TRUE)
  return(entropy)
}


#' Calculate the max entropy for the cells when given the same counts.
#'
#' @param x A sparse matrix which columns represent cells and rows represent features.
#'
#' @return A vector of the max cell entropy.
#'
#' @examples
#' x <- simCounts()
#' maxCE <- maxCellEntropy(x)
#' head(maxCE)
#' @importFrom Matrix t colSums
#' @export
maxCellEntropy <- function(x) {
  if (!class(x) %in% c("dgCMatrix", "dgTMatrix")) {
    stop("'x' must be sparse Matrix of class dgCMatrix or dgTMatrix")
  }
  nCount <- Matrix::colSums(x)
  nGene <- Matrix::colSums(x > 0)
  i <- nCount %/% nGene
  r <- nCount %% nGene
  pi <- i / nCount
  pr <- (i + 1) / nCount
  ei <- (nGene - r) * pi * log2(pi)
  er <- r * pr * log2(pr)
  ei[is.na(ei)] <- 0
  er[is.na(er)] <- 0
  entropy <- -(ei + er)
  return(entropy)
}

#' Calculate the entropy rate for the cells according to the feature counts.
#'
#' @param x A sparse matrix which columns represent cells and rows represent features.
#'
#' @return A vector of the cell entropy rate.
#'
#' @examples
#' x <- simCounts()
#' CER <- CellEntropyRate(x)
#' head(CER)
#' @export
CellEntropyRate <- function(x) {
  if (!class(x) %in% c("dgCMatrix", "dgTMatrix")) {
    stop("'x' must be sparse Matrix of class dgCMatrix or dgTMatrix")
  }
  er <- CellEntropy(x) / maxCellEntropy(x)
  return(er)
}
