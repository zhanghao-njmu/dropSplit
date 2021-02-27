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
#' @importFrom Matrix t colSums
#' @export
CellGini <- function(x, normalize = TRUE) {
  if (!class(x) %in% c("dgCMatrix", "dgTMatrix")) {
    stop("'x' must be sparse Matrix of class dgCMatrix or dgTMatrix")
  }
  p_matrix <- Matrix::t(Matrix::t(x) / Matrix::colSums(x))
  gini <- 1 - Matrix::colSums(p_matrix^2)
  if (isTRUE(normalize)) {
    nCount <- Matrix::colSums(x)
    nGene <- Matrix::colSums(x > 0)
    i <- nCount %/% nGene
    r <- nCount %% nGene
    pi <- i / nCount
    pr <- (i + 1) / nCount
    m <- (nGene - r) * pi^2 + r * pr^2
    gini <- gini + m
  }
  return(gini)
}

#' Calculate GiniScore under a threshold.
#'
#' @param x A vector of Gini index.
#' @param GiniThreshold A value used to estimate the feature-count redundancy from the Gini index. \code{x} larger than GiniThreshold will result in a GiniScore>0.5, else GiniScore<0.5.
#' @param group The group of \code{x}.
#' @return A vector of the GiniScore.
#'
#' @examples
#' x <- c(0.6, 0.7, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.98, 0.99)
#' GiniScore(x, 0.95)
#' @export
GiniScore <- function(x, GiniThreshold, group = NULL) {
  if (any(x < 0 | x > 1)) {
    stop("'x' must be in between 0 and 1.")
  }
  if (is.null(group)) {
    group <- rep(1, length(x))
  }
  r <- x - GiniThreshold
  score <- rep(0, length(x))
  for (g in unique(group)) {
    j <- group == g
    i <- r[j]
    score[j] <- ifelse(i >= 0, i / (max(x[j]) - GiniThreshold), i / (GiniThreshold - min(x[j])))
  }
  score[is.na(score)] <- 0
  score <- (score + 1) / 2
  return(score)
}
