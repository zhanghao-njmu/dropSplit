#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of three parts: quality control, model construction and droplet classification, summarizing features.
#' dropSplit provides some special QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent samples.
#' @param GiniThreshold A value used in \code{\link{GiniScore}} function. The higher, the more conservative dropSplit will be. Default is automatic.
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
#' @importFrom xgboost xgboost xgb.DMatrix xgb.importance xgb.dump
#' @importFrom Matrix t colSums
#' @importFrom methods as
#' @importFrom stats na.omit predict
#' @export
dropSplit <- function(counts, GiniThreshold = NULL, score_cutoff = 0.8, modelOpt = FALSE,
                      xgb_params = NULL, xgb_nrounds = 20, xgb_thread = 8,
                      bounds = list(),
                      xgb_nfold = 5, xgb_early_stopping_rounds = 3, xgb_metric = "auc",
                      opt_initPoints = length(bounds) + 1, opt_itersn = 10, opt_thread = 1, ...) {
  if (!class(counts) %in% c("matrix", "dgCMatrix", "dgTMatrix")) {
    stop("'counts' must be a dense matrix or a sparse matrix object.")
  }
  if (length(colnames(counts)) != ncol(counts) | length(rownames(counts)) != nrow(counts)) {
    stop("'counts' matrix must have both row(feature) names and column(cell) names.")
  }
  if (!is.null(GiniThreshold)) {
    if (GiniThreshold < 0 | GiniThreshold > 1) {
      stop("'GiniThreshold' must be between 0 and 1.")
    }
  }
  if (!is.null(score_cutoff)) {
    if (score_cutoff < 0 | score_cutoff > 1) {
      stop("'score_cutoff' must be between 0 and 1.")
    }
  }
  if (class(counts) == "matrix") {
    counts <- as(counts, "dgCMatrix")
  }

  message(">>> Start to define the credible cell-containing droplet...")
  set.seed(0)
  meta_info <- data.frame(row.names = colnames(counts))
  meta_info$nCount <- Matrix::colSums(counts)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  meta_info$nFeature <- Matrix::colSums(counts > 0)
  meta_info$nFeature_rank <- rank(-(meta_info$nFeature))

  if (min(meta_info$nCount) <= 0) {
    warning("'counts' has cells that nCount <=0. These cells will be remove in the following steps.",
      immediate. = TRUE
    )
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
  Cell_rank <- which.min(meta_info$RankMSE[1:inflection])
  Cell_count <- meta_info$nCount[Cell_rank]
  Uncertain_rank <- Cell_rank + which.max(meta_info$RankMSE[Cell_rank:length(meta_info$RankMSE)]) - 1
  Uncertain_count <- meta_info$nCount[Uncertain_rank]
  Empty_rank <- Uncertain_rank + find_peaks(-meta_info$RankMSE[Uncertain_rank:length(meta_info$RankMSE)],
    left_shoulder = round(0.01 * ncol(counts)),
    right_shoulder = round(0.05 * ncol(counts))
  )[1]
  Empty_count <- meta_info$nCount[Empty_rank]

  Cell_counts <- counts[, meta_info$nCount >= Cell_count]
  Uncertain_counts <- counts[, meta_info$nCount < Cell_count & meta_info$nCount >= Uncertain_count]
  Empty_counts <- counts[, meta_info$nCount < Uncertain_count & meta_info$nCount >= Empty_count]
  message(
    ">>> Summary of pre-defined droplets",
    "\n... Number of Cell: ", ncol(Cell_counts), "  Minimum nCounts: ", Cell_count,
    "\n... Number of Uncertain: ", ncol(Uncertain_counts), "  Minimum nCounts: ", Uncertain_count,
    "\n... Number of Empty: ", ncol(Empty_counts), "  Minimum nCounts: ", Empty_count,
    "\n... Number of Discarded: ", ncol(counts) - ncol(Cell_counts) - ncol(Uncertain_counts) - ncol(Empty_counts), "  Minimum nCounts: ", min(meta_info$nCount)
  )
  # ggplot(meta_info,aes(x=log10(1:nrow(meta_info)),y=RankMSE))+
  #   geom_point(alpha=0.1)+
  #   geom_vline(xintercept = log10(c(Uncertain_rank,Empty_rank)))

  message(">>> Simulate low depth cell-containing droplets from pre-defined Cell droplets...")
  i <- sample(x = 1:ncol(Cell_counts), size = ifelse(ncol(Empty_counts) < 100000, ncol(Empty_counts), 100000), replace = TRUE)
  Sim_counts <- Cell_counts[, i]
  colnames(Sim_counts) <- paste0("Sim-", 1:ncol(Sim_counts))
  Sim_nCount <- Matrix::colSums(Sim_counts)
  Sim_nCount_assign <- sample(unique(Matrix::colSums(Empty_counts)), ncol(Empty_counts), replace = TRUE)
  Sim_counts <- downsampleMatrix(x = Sim_counts, prop = Sim_nCount_assign / Sim_nCount, bycol = TRUE)

  message(">>> Calculate QC metrics for the droplets to be trained...")
  comb_counts <- cbind(Cell_counts, Sim_counts, Uncertain_counts, Empty_counts)
  comb_nCount <- Matrix::colSums(comb_counts)
  comb_nFeature <- Matrix::colSums(comb_counts > 0)
  comb_CellEntropy <- CellEntropy(comb_counts)
  comb_EntropyRate <- comb_CellEntropy / maxCellEntropy(comb_counts)
  comb_EntropyRate[is.na(comb_EntropyRate)] <- 1
  comb_Gini <- CellGini(comb_counts, normalize = T)
  if (is.null(GiniThreshold)) {
    GiniThreshold <- min(quantile(comb_Gini[colnames(Cell_counts)], 0.01), 0.99)
  }
  comb_GiniScore <- GiniScore(
    x = comb_Gini, GiniThreshold = GiniThreshold,
    group = c(rep("Cell", ncol(Cell_counts)), rep("Sim", ncol(Sim_counts)), rep("Uncertain", ncol(Uncertain_counts)), rep("Empty", ncol(Empty_counts)))
  )
  MTgene <- grep(x = rownames(comb_counts), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T, value = TRUE)
  RPgene <- grep(x = rownames(comb_counts), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T, value = TRUE)
  comb_MTprop <- Matrix::colSums(comb_counts[MTgene, ]) / Matrix::colSums(comb_counts)
  comb_RPprop <- Matrix::colSums(comb_counts[RPgene, ]) / Matrix::colSums(comb_counts)
  dat_other <- cbind(
    nCount = comb_nCount,
    nFeature = comb_nFeature,
    CellEntropy = comb_CellEntropy,
    EntropyRate = comb_EntropyRate,
    CellGini = comb_Gini,
    GiniScore = comb_GiniScore,
    MTprop = comb_MTprop,
    RPprop = comb_RPprop
  )
  rownames(dat_other) <- colnames(comb_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_other)
  ] <- dat_other[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]

  norm_counts <- Matrix::t(Matrix::t(comb_counts) / Matrix::colSums(comb_counts))
  dat <- cbind(Matrix::t(norm_counts), dat_other[, !colnames(dat_other) %in% c("CellGini", "GiniScore")])
  train <- dat[c(colnames(Cell_counts), colnames(Sim_counts), colnames(Empty_counts)), ]
  train_label <- c(rep(1, ncol(Cell_counts) + ncol(Sim_counts)), rep(0, ncol(Empty_counts)))
  to_predict <- dat[colnames(Uncertain_counts), ]

  if (isTRUE(modelOpt)) {
    opt <- xgbOptimization(
      dat = train, dat_label = train_label,
      bounds = bounds,
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
  message(">>> eXtreme Gradient Boosting(XGBoost) training with pre-defined classification...")
  xgb <- xgboost(
    data = xgb.DMatrix(data = train, label = train_label),
    nrounds = xgb_nrounds,
    early_stopping_rounds = xgb_early_stopping_rounds,
    params = xgb_params
  )
  if (nrow(to_predict) == 0) {
    message("No Uncertain droplets.")
    pred <- numeric(0)
  } else {
    pred <- predict(xgb, to_predict)
  }
  pred <- c(rep(1, ncol(Cell_counts)), pred)
  score <- pmin((pred + meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts)), "GiniScore"]) / 2, pred)

  meta_info[, "preDefinedClass"] <- "Discarded"
  meta_info[colnames(Cell_counts), "preDefinedClass"] <- "Cell"
  meta_info[colnames(Uncertain_counts), "preDefinedClass"] <- "Uncertain"
  meta_info[colnames(Empty_counts), "preDefinedClass"] <- "Empty"
  meta_info[, "dropSplitClass"] <- meta_info[, "preDefinedClass"]
  meta_info[, "dropSplitScore"] <- 0
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts)), "dropSplitScore"] <- score
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts)), "dropSplitClass"] <- ifelse(
    score > score_cutoff, "Cell", ifelse(score < 1 - score_cutoff, "Empty", "Uncertain")
  )

  meta_info <- meta_info[raw_cell_order, ]
  message(
    ">>> Summary of dropSplit-defined droplet:",
    "\n... Number of Cell: ", sum(meta_info$dropSplitClass == "Cell"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Cell", "nCount"]),
    "\n... Number of Uncertain: ", sum(meta_info$dropSplitClass == "Uncertain"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Uncertain", "nCount"]),
    "\n... Number of Empty: ", sum(meta_info$dropSplitClass == "Empty"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Empty", "nCount"]),
    "\n... Number of Discarded: ", sum(meta_info$dropSplitClass == "Discarded"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Discarded") > 1, min(meta_info[meta_info$dropSplitClass == "Discarded", "nCount"]), min(meta_info[, "nCount"])),
    "\n>>> Pre-defined as 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell"),
    "\n... Mean CellGini:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "CellGini"]), 3),
    "\n... Mean dropSplitScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "dropSplitScore"]), 3),
    "\n>>> Pre-defined as 'Uncertain' switch to 'Cell': ", sum(meta_info$preDefinedClass == "Uncertain" & meta_info$dropSplitClass == "Cell"),
    "\n... Mean CellGini:", round(mean(meta_info[meta_info$preDefinedClass == "Uncertain" & meta_info$dropSplitClass == "Cell", "CellGini"]), 3),
    "\n... Mean dropSplitScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3)
  )
  importance_matrix <- xgb.importance(model = xgb)
  tree <- xgb.dump(xgb, with_stats = TRUE)
  result <- list(
    meta_info = meta_info,
    train = train, train_label = train_label,
    model = xgb,
    importance_matrix = importance_matrix,
    tree = tree
  )

  # xgb.ggplot.importance(importance_matrix[1:20, ])
  # meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  # meta_info2 <- meta_info[meta_info$dropSplitClass=="Cell",]
  # meta_info2$nCount_rank <- rank(-(meta_info2$nCount))
  # meta_info2$nFeature_rank <- rank(-(meta_info2$nFeature))
  # meta_info2 <- meta_info2[order(meta_info2$nCount_rank, decreasing = FALSE), ]
  # meta_info2$RankSE <- (log10(meta_info2$nCount_rank) - log10(meta_info2$nFeature_rank))^2
  # x <- meta_info2$RankSE
  # for (t in seq_len(2)) {
  #   x <- runMean(x, n = 100)
  #   x[is.na(x)] <- na.omit(x)[1]
  # }
  # meta_info2$RankMSE <- x
  # ggplot(meta_info,aes(x=log10(1:nrow(meta_info)),y=RankMSE,color=dropSplitClass))+geom_point(alpha=0.1)+
  #   geom_point(data = meta_info2,aes(x=log10(1:nrow(meta_info2)),y=RankMSE,color=dropSplitClass),alpha=0.1,color="yellow",shape=21)
  # ggplot(meta_info,aes(x=log10(nCount),y=CellEntropy,color=dropSplitClass))+
  #   geom_point(alpha=0.1)
  # geom_point(data = meta_info[meta_info$CellGini>0.99,])+

  return(result)
}
