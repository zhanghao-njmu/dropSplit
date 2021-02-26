#' Bayesian optimization for XGBoost.
#'
#' @description Maximizes a xgboost evaluation metric within a set of bounds. After the function is sampled a pre-determined number of times, a Gaussian process is fit to the results. An acquisition function is then maximized to determine the most likely location of the global maximum of the user defined XGBoost evaluation metric. This process is repeated for a set number of iterations.
#' @param dat A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent samples.
#' @param dat_label A vector of response classification values.
#' @param bounds A named list of lower and upper bounds for \code{params} in \code{\link[xgboost]{xgb.cv}}. The names of the list should be arguments passed to xgb.cv Use "L" suffix to indicate integers. A fixed parameter should be a two-length vector with the same value, i.e. bound=list(lambda = c(5, 5))
#' @param xgb_nfold The original dataset is randomly partitioned into nfold equal size subsamples.
#' @param xgb_nround Max number of boosting iterations.
#' @param xgb_early_stopping_rounds If NULL, the early stopping function is not triggered. If set to an integer k, training with a validation set will stop if the performance doesn't improve for k rounds. Setting this parameter engages the \code{cb.early.stop} callback.
#' @param xgb_metric A evaluation metric to be used in cross validation and will to be maximized. Possible options are:
#' \itemize{
#' \item \code{auc} Area under curve
#' \item \code{aucpr} Area under PR curve
#' }
#' @param xgb_thread Number of thread used in \code{\link[xgboost]{xgb.cv}}.
#' @param opt_initPoints Number of points to initialize the process with. Points are chosen with latin hypercube sampling within the bounds supplied.
#' @param opt_itersn The total number of times \code{xgb.cv} will be run after initialization.
#' @param opt_thread Number of thread used in \code{\link[ParBayesianOptimization]{bayesOpt}}.
#' @param ... Other arguments passed to \code{\link[ParBayesianOptimization]{bayesOpt}}.
#'
#' @return A list of two object:
#' \describe{
#' \item{bayesOpt}{An object of class bayesOpt containing information about the process.}
#' \item{BestPars}{A list containing the parameters which resulted in the highest returned Score.}
#' }
#'
#' @examples
#' library("xgboost")
#' data(agaricus.train, package = "xgboost")
#' dat <- agaricus.train$data
#' dat_label <- agaricus.train$label
#' bounds <- list(max_depth = c(1L, 5L), min_child_weight = c(0, 25), subsample = c(0.25, 1))
#' result <- xgbOptimization(dat = dat, dat_label = dat_label, bounds = bounds, opt_thread = 2)
#' result
#' @importFrom ParBayesianOptimization bayesOpt getBestPars
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom BiocGenerics clusterExport clusterEvalQ
#' @importFrom xgboost xgb.cv
#' @export
xgbOptimization <- function(dat, dat_label, bounds = list(),
                            xgb_nfold = 5, xgb_nround = 20, xgb_early_stopping_rounds = 5, xgb_metric = "auc", xgb_thread = 8,
                            opt_initPoints = length(bounds) + 1, opt_itersn = 10, opt_thread = 1, ...) {
  set.seed(0)
  if (length(bounds) == 0) {
    bounds <- list(
      eta = c(0.01, 0.3),
      gamma = c(0L, 10L),
      max_depth = c(10L, 30L),
      min_child_weight = c(1, 10),
      subsample = c(0.5, 1),
      max_delta_step = c(0.7, 1),
      alpha = c(1, 10),
      lambda = c(1, 10)
    )
  }
  xgbScoreFun <- defineScoreFun(
    dat = dat, dat_label = dat_label,
    nfold = xgb_nfold, nround = xgb_nround, early_stopping_rounds = xgb_early_stopping_rounds,
    metric = xgb_metric, xgb_thread = xgb_thread
  )

  if (opt_thread > 1) {
    cl <- makeCluster(opt_thread)
    registerDoParallel(cl)
    clusterExport(cl, c(
      paste0(substitute(dat)),
      paste0(substitute(dat_label))
    ))
    clusterEvalQ(cl, expr = {
      library(xgboost)
    })

    optObj <- bayesOpt(
      FUN = xgbScoreFun,
      bounds = bounds,
      initPoints = opt_initPoints,
      iters.n = opt_itersn,
      iters.k = opt_thread,
      parallel = TRUE,
      ...
    )

    stopCluster(cl)
    registerDoSEQ()
  } else {
    optObj <- bayesOpt(
      FUN = xgbScoreFun,
      bounds = bounds,
      initPoints = opt_initPoints,
      iters.n = opt_itersn,
      iters.k = 1,
      parallel = FALSE,
      ...
    )
  }

  result <- list(bayesOpt = optObj, BestPars = getBestPars(optObj))
  return(result)
}

defineScoreFun <- function(dat, dat_label,
                           ...,
                           nfold = 5, nround = 20, early_stopping_rounds = 5, metric = "auc",
                           xgb_thread = 8) {
  xgbScore <- function(.dat = get("dat"), .dat_label = get("dat_label"),
                       ...,
                       .xgb_nfold = get("xgb_nfold"), .xgb_nround = get("xgb_nround"), .xgb_early_stopping_rounds = get("xgb_early_stopping_rounds"), .xgb_metric = get("xgb_metric"),
                       .xgb_thread = get("xgb_thread")) {
    xgbDMatrix <- xgb.DMatrix(data = .dat, label = .dat_label)
    Pars <- list(
      booster = "gbtree",
      ...,
      objective = "binary:logistic",
      nthread = .xgb_thread
    )
    xgbcv <- xgb.cv(
      params = Pars,
      data = xgbDMatrix,
      nround = .xgb_nround,
      nfold = .xgb_nfold,
      early_stopping_rounds = .xgb_early_stopping_rounds,
      metrics = list(.xgb_metric),
      prediction = TRUE,
      showsd = TRUE,
      maximize = TRUE,
      verbose = 0
    )
    return(
      list(
        Score = max(xgbcv$evaluation_log[[paste0("test_", metric, "_mean")]]),
        nrounds = xgbcv$best_iteration
      )
    )
  }
  return(xgbScore)
}
