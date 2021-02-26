#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of three parts: quality control, model construction and droplet classification, summarizing features.
#' dropSplit also provides some special QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also implements a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent samples.
#' @param score_cutoff A cutoff value for dropSplit to determine if a droplet is cell-containing or empty.
#' @param modelOpt Whether to optimize the model using \code{\link{xgbOptimization}}.If true, will overwrite the parameters list in \code{xgb_params}.
#' @param xgb_params The \code{list} of XGBoost parameters.
#' @inheritParams xgbOptimization
#' @param ... Other arguments passed to \code{\link{xgbOptimization}}.
#'
#' @return A list of two object:
#' \describe{
#' \item{meta_info}{A data.frame object of evaluation metrics to be used in dropSplit and the final droplet classification.}
#' \item{model}{The XGBoost model used in dropSplit for classification.}
#' }
#'
#' @examples
#' counts <- DropletUtils:::simCounts()
#' colnames(counts) <- paste0("Cell-", 1:ncol(counts))
#' result <- dropSplit(counts)
#' head(result$meta_info)
#' @importFrom TTR runMean
#' @importFrom DropletUtils downsampleMatrix
#' @importFrom xgboost xgboost xgb.DMatrix
#' @importFrom Matrix t colSums
#' @importFrom methods as
#' @importFrom stats na.omit predict
#' @export
dropSplit <- function(counts, score_cutoff = 0.8, modelOpt = FALSE,
                      xgb_params = NULL, xgb_nrounds = 20, xgb_thread = 8,
                      bounds = list(),
                      xgb_nfold = 5, xgb_early_stopping_rounds = 3, xgb_metric = "auc",
                      opt_initPoints = length(bounds) + 1, opt_itersn = 10, opt_thread = 1, ...) {
  if (!class(counts) %in% c("matrix", "dgCMatrix", "dgTMatrix")) {
    stop("'counts' must be a dense matrix or a sparse matrix object.")
  }
  if (class(counts) == "matrix") {
    counts <- as(counts, "dgCMatrix")
  }
  if (length(colnames(counts)) != ncol(counts) | length(rownames(counts)) != nrow(counts)) {
    stop("'counts' matrix must have both row(feature) names and column(cell) names")
  }

  message(">>> Start to define the credible cell-containing droplet...\n")
  set.seed(0)
  meta_info <- data.frame(row.names = colnames(counts))
  meta_info$nCount <- Matrix::colSums(counts)
  meta_info$nCount_rank <- rank(-(Matrix::colSums(counts)))
  meta_info$nFeature <- Matrix::colSums(counts > 0)
  meta_info$nFeature_rank <- rank(-(Matrix::colSums(counts > 0)))

  if (min(meta_info$nCount) <= 0) {
    warning("'counts' has cells that nCount <=0. These cells will be remove in the following steps.\n",
            immediate. = TRUE,noBreaks. = TRUE)
    meta_info <- meta_info[meta_info$nCount > 0, ]
  }
  raw_cell_order <- rownames(meta_info)
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  counts <- counts[, rownames(meta_info)]

  nCount_inflection <- find_inflection(meta_info$nCount)$index[1]
  nFeature_inflection <- find_inflection(meta_info$nFeature)$index[1]
  inflection <- min(nCount_inflection, nFeature_inflection)
  meta_info$RankSE <- (log10(meta_info$nCount_rank) - log10(meta_info$nFeature_rank))^2
  x <- meta_info$RankSE
  for (t in seq_len(2)) {
    x <- runMean(x, n = 100)
    x[is.na(x)] <- na.omit(x)[1]
  }
  meta_info$RankMSE <- x
  certain_rank <- which.min(meta_info$RankMSE[1:inflection])
  certain_count <- meta_info$nCount[certain_rank]
  uncertain_rank <- inflection + find_peaks(-meta_info$RankMSE[inflection:length(meta_info$RankMSE)])[1]
  uncertain_count <- meta_info$nCount[uncertain_rank]

  certain_cells <- counts[, meta_info$nCount >= certain_count]
  certain_nCount <- Matrix::colSums(certain_cells)
  uncertain_cells <- counts[, meta_info$nCount < certain_count & meta_info$nCount >= uncertain_count]
  uncertain_nCount <- Matrix::colSums(uncertain_cells)
  message(
    ">>> Summary of pre-defined droplet",
    "\n... Number of Cell: ", ncol(certain_cells), "  Minimum nCounts: ", certain_count,
    "\n... Number of Empty: ", ncol(uncertain_cells), "  Minimum nCounts: ", uncertain_count,
    "\n... Number of Discarded: ", ncol(counts) - ncol(certain_cells) - ncol(uncertain_cells), "  Minimum nCounts: ", min(meta_info$nCount), "\n"
  )

  message(">>> Simulate the low depth cell-containing droplet from the pre-defined Cell...\n")
  i <- sample(x = 1:ncol(certain_cells), size = ncol(uncertain_cells), replace = TRUE)
  sim_cells <- certain_cells[, i]
  colnames(sim_cells) <- paste0("sim_cells-", 1:ncol(sim_cells))
  sim_nCount <- Matrix::colSums(sim_cells)
  sim_nCount_assign <- sample(unique(uncertain_nCount), ncol(uncertain_cells), replace = TRUE)
  sim_cells <- downsampleMatrix(x = sim_cells, prop = sim_nCount_assign / sim_nCount, bycol = TRUE)

  comb_cells <- cbind(certain_cells, sim_cells, uncertain_cells)
  dat <- Matrix::t(comb_cells)
  dat_label <- c(rep(1, ncol(certain_cells) + ncol(sim_cells)), rep(0, ncol(uncertain_cells)))

  message(">>> Calculate QC metrics for the droplets to be trained...\n")
  comb_CellEntropy <- CellEntropy(comb_cells)
  comb_EntropyRate <- comb_CellEntropy / maxCellEntropy(comb_cells)
  comb_EntropyRate[is.na(comb_EntropyRate)] <- 1
  comb_gini <- CellGini(comb_cells, normalize = T)
  comb_nCount <- Matrix::colSums(comb_cells)
  comb_nFeature <- Matrix::colSums(comb_cells > 0)
  MTgene <- grep(x = rownames(comb_cells), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T, value = TRUE)
  RPgene <- grep(x = rownames(comb_cells), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T, value = TRUE)
  comb_MTprop <- Matrix::colSums(comb_cells[MTgene, ]) / Matrix::colSums(comb_cells)
  comb_RPprop <- Matrix::colSums(comb_cells[RPgene, ]) / Matrix::colSums(comb_cells)
  dat_other <- cbind(
    CellEntropy = comb_CellEntropy,
    EntropyRate = comb_EntropyRate,
    CellGini = comb_gini,
    nCount = comb_nCount,
    nFeature = comb_nFeature,
    MTprop = comb_MTprop,
    RPprop = comb_RPprop
  )
  rownames(dat_other) <- colnames(comb_cells)
  meta_info <- merge(
    x = meta_info, y = dat_other[, c("CellEntropy", "EntropyRate", "CellGini", "MTprop", "RPprop")],
    by = "row.names", all.x = TRUE, all.y = FALSE,
  )
  rownames(meta_info) <- meta_info$Row.names
  meta_info <- meta_info[, -1]
  dat <- cbind(dat, dat_other)

  if (isTRUE(modelOpt)) {
    opt <- xgbOptimization(dat, dat_label,
      bounds = list(),
      xgb_nfold = xgb_nfold, xgb_nrounds = xgb_nrounds, xgb_early_stopping_rounds = xgb_early_stopping_rounds, xgb_metric = xgb_metric, xgb_thread = xgb_thread,
      opt_initPoints = opt_initPoints, opt_itersn = opt_itersn, opt_thread = opt_thread, ...
    )
    xgb_params <- c(opt$BestPars,
      eval_metric = "error",
      eval_metric = "auc",
      objective = "binary:logistic",
      nthread = xgb_thread
    )
  }
  if (is.null(xgb_params)) {
    xgb_params <- list(
      eta = 0.3,
      gamma = 1,
      max_depth = 20,
      min_child_weight = 7,
      subsample = 0.7,
      max_delta_step = 1,
      alpha = 5,
      lambda = 10,
      eval_metric = "error",
      eval_metric = "auc",
      objective = "binary:logistic",
      nthread = xgb_thread
    )
  }
  message(">>> eXtreme Gradient Boosting(XGBoost) training for the pre-defined the droplets...\n")
  xgb <- xgboost(
    data = xgb.DMatrix(data = dat, label = dat_label),
    nrounds = xgb_nrounds,
    early_stopping_rounds = xgb_early_stopping_rounds,
    params = xgb_params
  )
  pred <- predict(xgb, dat[rownames(meta_info), ])
  meta_info[, "preDefinedClass"] <- "Discarded"
  meta_info[colnames(certain_cells), "preDefinedClass"] <- "Cell"
  meta_info[colnames(uncertain_cells), "preDefinedClass"] <- "Empty"
  meta_info[, "dropSplitClass"] <- "Discarded"
  meta_info[, "dropSplitScore"] <- -1
  meta_info[rownames(meta_info), "dropSplitScore"] <- pred
  meta_info[rownames(meta_info), "dropSplitClass"] <- ifelse(pred > score_cutoff, "Cell", "Empty")
  meta_info <- meta_info[raw_cell_order, ]
  message(
    ">>> Summary of dropSplit-defined droplet",
    "\n... Number of Cell: ", sum(meta_info$dropSplitClass == "Cell"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Cell", "nCount"]),
    "\n... Number of Empty: ", sum(meta_info$dropSplitClass == "Empty"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Empty", "nCount"]),
    "\n... Number of Discarded: ", sum(meta_info$dropSplitClass == "Discarded"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Discarded") > 1, min(meta_info[meta_info$dropSplitClass == "Discarded", "nCount"]), min(meta_info[, "nCount"])), "\n"
  )
  result <- list(meta_info = meta_info, model = xgb)
  return(result)
}
