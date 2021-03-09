#' Find the index of peaks in a sequence of numeric values.
#' @description Find the index of peaks in a sequence of numeric values.
#' Peaks value should be larger than any values among the left/right shoulder.
#' If no peaks found, function will return the index of the max number.
#'
#' @param x A numeric vector.
#' @param left_shoulder A integer value. Peaks value should larger than any value in the left \code{left_shoulder} value.
#' @param right_shoulder A integer value. Peaks value should larger than any value in the right \code{right_shoulder} value.
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
#' x <- simSimpleCounts()
#' inflection <- find_inflection(Matrix::colSums(x))
#' inflection
#' @importFrom stats smooth.spline predict
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
  curvature <- curvatureCalcluate(fitted_y, fitted_x)$curvature
  inflection_y <- fitted_y[which.max(curvature)]
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

#' Calculate Score under a threshold within groups.
#'
#' @param x A vector of metric used to score.
#' @param threshold A value used to calculate the score. \code{x} larger than threshold will result in a Score>0.5, else Score<0.5.
#' @param group The groups of \code{x}. If \code{NULL}, elements in \code{x} will be treated as the same group.
#' @param upper,lower A named vector of value to specify the upper/lower value for each group. Default is to find the upper for each group automatically.
#' @param higher_score,lower_score If provided, function will return a fixed score. \code{x>threshold} return \code{higher_score}, else return \code{higher_score}.
#' @return A vector of the Score. A Score>0.5 represent that corresponding x is larger than the \code{threshold} when fuzz=\code{FALSE}.
#'
#' @examples
#' x <- c(0.6, 0.7, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.98, 0.99)
#' Score(x, 0.95)
#' @importFrom stats quantile
#' @export
Score <- function(x, threshold, group = NULL, upper = NULL, lower = NULL, higher_score = NULL, lower_score = NULL) {
  if (is.null(group)) {
    group <- rep(1, length(x))
  }
  if (!is.null(upper)) {
    if (all(upper < threshold)) {
      stop("values in the 'upper' must be >=threshold")
    }
  }
  if (!is.null(lower)) {
    if (all(lower > threshold)) {
      stop("values in the 'lower' must be <=threshold")
    }
  }
  r <- x - threshold
  score <- rep(0, length(x))
  for (g in unique(group)) {
    j <- which(group == g)
    i <- r[j]

    ilarge <- i[which(i >= 0)]
    if (length(ilarge) > 0) {
      if (!g %in% (names(upper))) {
        maxi <- quantile(ilarge, 0.9)
      } else {
        maxi <- upper[[g]] - threshold
      }
      score[j[which(i >= 0)]] <- ifelse(maxi == 0 & ilarge != 0, 1, ifelse(maxi == 0 & ilarge == 0, 0, ilarge / maxi))
    }

    ilow <- i[which(i < 0)]
    if (length(ilow) > 0) {
      if (!g %in% (names(lower))) {
        mini <- quantile(ilow, 0.1)
      } else {
        mini <- threshold - lower[[g]]
      }
      score[j][which(i < 0)] <- ifelse(mini == 0 & ilow != 0, -1, ifelse(mini == 0 & ilow == 0, 0, -ilow / mini))
    }
  }
  score <- (score + 1) / 2
  score[score > 1] <- 1
  score[score < 0] <- 0
  if (!is.null(higher_score)) {
    score[score > 0.5] <- higher_score
  }
  if (!is.null(lower_score)) {
    score[score <= 0.5] <- lower_score
  }
  return(score)
}

