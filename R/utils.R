#' Find the index of peaks in a sequence of numeric values.
#' @description Find the index of peaks in a sequence of numeric values.
#' Peaks value should be larger than any values among the left/right shoulder.
#' If no peaks found, function will return the index of the max number.
#'
#' @param x A numeric vector.
#' @param left_shoulder A integer value.Peaks value should larger than any value in the left \code{left_shoulder} value.
#' @param right_shoulder A integer value.Peaks value should larger than any value in the right \code{right_shoulder} value.
#'
#' @return A vector of the index of peaks. If no peaks found,return the index of the max number.
#'
#' @examples
#' x <- c(0, 1, 3, 6, 9, 12, 11, 7, 9, 5, 1, 9, 0, 1, 2)
#' pks <- find_peaks(x, left_shoulder = 3, right_shoulder = 3)
#' pks
find_peaks <- function(x, left_shoulder = 5000, right_shoulder = 20000) {
  raw_order <- 1:length(x)
  na_order <- which(is.na(x))
  if (length(na_order) == 0) {
    filter_order <- raw_order
  } else {
    x <- na.omit(x)
    filter_order <- raw_order[-na_order]
  }
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i) {
    f <- i - left_shoulder + 1
    b <- i + right_shoulder + 1
    if (f <= 0 | b > length(x)) {
      return(numeric(0))
    }
    if (all(x[c(f:i, (i + 2):b)] <= x[i + 1])) {
      return(i + 1)
    } else {
      return(numeric(0))
    }
  })
  pks_order <- filter_order[unlist(pks)]
  if (length(pks_order) == 0) {
    pks_order <- which.max(x)
  }
  return(pks_order)
}

#' Find the index of inflection in a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return A list of the index of the inflection and the corresponding value.
#'
#' @examples
#' x <- DropletUtils:::simCounts()
#' inflection <- find_inflection(Matrix::colSums(x))
#' inflection
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
  if (min(curvature) < 0 & abs(min(curvature)) > mean(curvature) + 2 * sd(curvature)) {
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
curvatureCalcluate <- function(x, y) {
  d1n <- diff(y) / diff(x)
  d2n <- diff(d1n) / diff(x[-1] - diff(x) / 2)
  d1n <- c(d1n, d1n[length(d1n)])
  d2n <- c(d2n[1], d2n, d2n[length(d2n)])
  curvature <- d2n / (1 + d1n^2)^1.5
  return(list(x = x, y = y, d1n = d1n, d2n = d2n, curvature = curvature))
}
