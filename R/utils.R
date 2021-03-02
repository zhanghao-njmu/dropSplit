#' Find the index of peaks in a sequence of numeric values.
#' @description Find the index of peaks in a sequence of numeric values.
#' Peaks value should be larger than any values among the left/right shoulder.
#' If no peaks found, function will return the index of the max number.
#'
#' @param x A numeric vector.
#' @param left_shoulder A integer value.Peaks value should larger than any value in the left \code{left_shoulder} value.
#' @param right_shoulder A integer value.Peaks value should larger than any value in the right \code{right_shoulder} value.
#'
#' @return A vector of the index of peaks.
#'
#' @examples
#' x <- c(0, 1, 3, 6, 9, 12, 11, 7, 9, 5, 1, 9, 0, 1, 2)
#' pks <- find_peaks(x, left_shoulder = 3, right_shoulder = 3)
#' pks
#' @importFrom stats na.omit
#' @export
find_peaks <- function(x, left_shoulder = 10000, right_shoulder = 10000) {
  raw_order <- 1:length(x)
  na_order <- which(is.na(x))
  if (length(na_order) == 0) {
    filter_order <- raw_order
  } else {
    x <- na.omit(x)
    filter_order <- raw_order[-na_order]
  }
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  shape <- c(-2, shape, -2)
  pks <- sapply(which(shape < 0), FUN = function(i) {
    f <- i - left_shoulder
    b <- i + right_shoulder
    f <- max(1, f)
    b <- min(length(x), b)
    pk <- i
    pkleft <- max(1, i - 1)
    pkright <- min(length(x), i + 1)
    # if (f <= 0 | b > length(x)) {
    #   return(numeric(0))
    # }
    if (all(x[c(f:pkleft, pkright:b)] <= x[pk]) & any(x[c(f:pkleft, pkright:b)] < x[pk])) {
      return(pk)
    } else {
      return(numeric(0))
    }
  })
  pks_index <- filter_order[unlist(pks)]
  return(pks_index)
}

#' Find the index of inflection in a numeric vector.
#'
#' @param x A numeric vector.
#' @param df The desired equivalent number of degrees of freedom (trace of the smoother matrix).
#'
#' @return A list of the index of the inflection and the corresponding value.
#'
#' @examples
#' x <- simCounts()
#' inflection <- find_inflection(Matrix::colSums(x))
#' inflection
#' @importFrom stats smooth.spline predict sd median quantile
#' @export
find_inflection <- function(x, df = 20) {
  r <- x
  x <- x[x > 0]
  x <- sort(x, decreasing = TRUE)
  raw_x <- rank(-x)
  raw_y <- x
  sp_fit1 <- smooth.spline(x = log10(raw_x), y = log10(raw_y), df = df)
  fitted <- predict(sp_fit1)
  fitted_x <- 10^fitted$x
  fitted_y <- 10^fitted$y
  curvature <- curvatureCalcluate(fitted_x, fitted_y)$curvature
  if (min(curvature) < 0 & min(curvature) < median(curvature) - 1.5 * (quantile(curvature, 0.75) - quantile(curvature, 0.25))) {
    inflection_y <- fitted_y[which(curvature == max(curvature[which.min(curvature):length(curvature)]))]
  } else {
    inflection_y <- fitted_y[which.max(curvature)]
  }
  value <- max(r[r < inflection_y])
  index <- which(r == value)
  return(list(index = index, value = value))
}

#' Calculate curvature from the given x and y.
#'
#' @param x,y A numeric vector.
#'
#' @return A list of the index of the raw x, y, and the corresponding first derivative, second derivative and curvature.
#'
#' @examples
#' y <- c(0, 1, 3, 6, 9, 12, 11, 7, 9, 5, 1, 9, 0, 1, 2)
#' x <- 1:length(y)
#' result <- curvatureCalcluate(x = x, y = y)
#' result
#' @export
curvatureCalcluate <- function(x, y) {
  d1n <- diff(y) / diff(x)
  d2n <- diff(d1n) / diff(x[-1] - diff(x) / 2)
  d1n <- c(d1n, d1n[length(d1n)])
  d2n <- c(d2n[1], d2n, d2n[length(d2n)])
  curvature <- d2n / (1 + d1n^2)^1.5
  return(list(x = x, y = y, d1n = d1n, d2n = d2n, curvature = curvature))
}

#' Simulate counts for the empty and cell-containing droplets.
#'
#' @param ngenes Total gene number for all the simulated counts.
#' @param nempty,nlarge,nsmall Empty, large cell and small cell droplets number.
#' @param empty.prof,large.prof,small.prof The overall gene expression profile distribution in distinct type of droplets.
#' @param empty.rate rate parameters in the \code{\link[stats]{Exponential}} function for empty droplets simulation.
#' @param large.rate,small.rate rate parameters in the \code{\link[stats]{GammaDist}} for large or small cell simulation.
#' @param large.shape,small.shape shape parameters in the \code{\link[stats]{GammaDist}} for large or small cell simulation.
#' @param RemoveZeroCol Whether to remove all zero-valued columns.
#'
#' @return A sparse Matrix of class "dgCMatrix".
#'
#' @examples
#' counts <- simCounts()
#' counts
#' @importFrom Matrix colSums
#' @importFrom stats rexp rgamma rpois runif
#' @export
simCounts <- function(ngenes = 5000, nempty = 20000, nlarge = 2000, nsmall = 200,
                      empty.prof = seq_len(ngenes), empty.rate = 0.04,
                      large.prof = empty.prof, large.rate = 0.01, large.shape = 10,
                      small.prof = runif(ngenes), small.rate = 0.1, small.shape = 20,
                      RemoveZeroCol = TRUE) {
  empty.prof <- empty.prof / sum(empty.prof)
  large.prof <- large.prof / sum(large.prof)
  small.prof <- small.prof / sum(small.prof)
  total.count <- rexp(nempty, rate = empty.rate)
  empty.counts <- matrix(rpois(ngenes * nempty, lambda = outer(
    empty.prof,
    total.count
  )), ncol = nempty, nrow = ngenes)
  empty.counts <- as(empty.counts, "dgCMatrix")
  total.count <- rgamma(nlarge, shape = large.shape, rate = large.rate)
  large.counts <- matrix(rpois(ngenes * nlarge, lambda = outer(
    large.prof,
    total.count
  )), ncol = nlarge, nrow = ngenes)
  large.counts <- as(large.counts, "dgCMatrix")
  total.count <- rgamma(nsmall, shape = small.shape, rate = small.rate)
  small.counts <- matrix(rpois(ngenes * nsmall, lambda = outer(
    small.prof,
    total.count
  )), ncol = nsmall, nrow = ngenes)
  small.counts <- as(small.counts, "dgCMatrix")
  out <- cbind(empty.counts, large.counts, small.counts)
  colnames(out) <- c(
    paste0("Empty-", seq_len(nempty)),
    paste0("LargeCell-", seq_len(nlarge)),
    paste0("SmallCell-", seq_len(nsmall))
  )
  rownames(out) <- paste0("Gene-", seq_len(ngenes))
  if (isTRUE(RemoveZeroCol)) {
    out <- out[, Matrix::colSums(out) > 0]
  }
  return(out)
}
