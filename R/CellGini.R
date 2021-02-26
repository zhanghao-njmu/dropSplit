#' Calculate the Gini index for the cells according to the feature counts.
#'
#' @param x A sparse matrix which columns represent cells and rows represent features.
#' @param normalize If \code{TRUE}, CellGini will add a number which represent the max Gini value when given the same counts.
#'
#' @return A vector of the cell entropy rate.
#'
#' @examples
#' x <- DropletUtils:::simCounts()
#' CG <- CellGini(x, normalize = TRUE)
#' head(CG)
#' @importFrom Matrix colSums
CellGini <- function(x, normalize = TRUE) {
  if (!class(x) %in% c("dgCMatrix", "dgTMatrix")) {
    stop("'x' must be sparse Matrix of class dgCMatrix or dgTMatrix")
  }
  p_matrix <- t(t(x) / colSums(x))
  gini <- 1 - colSums(p_matrix^2)
  if (isTRUE(normalize)) {
    nCount <- colSums(x)
    nGene <- colSums(x > 0)
    i <- nCount %/% nGene
    r <- nCount %% nGene
    pi <- i / nCount
    pr <- (i + 1) / nCount
    m <- (nGene - r) * pi^2 + r * pr^2
    gini <- gini + m
  }
  return(gini)
}
