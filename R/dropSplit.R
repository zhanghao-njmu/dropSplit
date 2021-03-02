#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of three parts: quality control, model construction and droplet classification, summarizing features.
#' dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent samples.
#' @param score_cutoff A cutoff value of dropSplitScore to determine if a droplet is cell-containing or empty. Default is 0.9, i.e. a droplet with dropSplitScore>0.9 will be classified as \code{Cell} and a with dropSplitScore<0.1  will be classified as \code{Empty}.
#' @param GiniThreshold A value used in \code{\link{GiniScore}} function. The higher, the more conservative dropSplit will be. Default is automatic.
#' @param Uncertain_downsample Whether use a downsampled Uncertain droplets for predicting. Default is FALSE.
#' @param Uncertain_downsample_times Number of downsample times for each Uncertain droplet. dropSplitScore for downsampled droplets from the same Uncertain droplet will be averaged. Default is 6.
#' @param predict_Uncertain_only Whether to predict only the Uncertain droplets. Default is \code{TRUE}.
#' @param Cell_rank,Uncertain_rank,Empty_rank Custom nCount_Rank value to mark the droplets as Cell, Uncertain and Empty labels for the data to be trained. Default is automatic. But useful when the default value is considered to be wrong from the RankMSE plot.
#' @param modelOpt Whether to optimize the model using \code{\link{xgbOptimization}}.If true, will overwrite the parameters list in \code{xgb_params}.
#' @param xgb_params The \code{list} of XGBoost parameters.
#' @inheritParams xgbOptimization
#' @param ... Other arguments passed to \code{\link{xgbOptimization}}.
#'
#' @return A list of seven objects:
#' \describe{
#' \item{meta_info}{A \code{data.frame} object of evaluation metrics to be used in dropSplit and the final droplet classification.}
#' \item{train}{The dataset trained in the XGBoost model. It consists of two pre-defined droplets: Cell(real-world + simulated) and Empty.}
#' \item{train_label}{Labels for the \code{train}. 0 represents 'Empty', 1 represents 'Cell'.}
#' \item{to_predict}{The dataset that to be predicted. It consists of all three pre-defined droplets: Cell, Uncertain and Empty.}
#' \item{model}{The XGBoost model used in dropSplit for classification.}
#' \item{importance_matrix}{A \code{data.frame} of feature importances in the classification model.}
#' \item{tree}{The trees from the classification model.}
#' }
#'
#' @examples
#' # Simulate a counts matrix including 20000 empty droplets, 2000 large cells and 200 small cells.
#' counts <- simCounts(nempty = 20000, nlarge = 2000, nsmall = 200)
#' counts_label <- gsub(pattern = "-.*", replacement = "", x = colnames(counts), perl = TRUE)
#' result <- dropSplit(counts)
#' head(result$meta_info)
#'
#' dropSplitClass <- result$meta_info$dropSplitClass
#' # True positive
#' sum(counts_label %in% c("LargeCell", "SmallCell") & dropSplitClass == "Cell")
#' # False negative
#' sum(counts_label %in% c("LargeCell", "SmallCell") & dropSplitClass != "Cell")
#' # True negative
#' sum(counts_label == "Empty" & dropSplitClass != "Cell")
#' # False positive
#' sum(counts_label == "Empty" & dropSplitClass == "Cell")
#' @importFrom TTR runMean
#' @importFrom DropletUtils downsampleMatrix
#' @importFrom xgboost xgboost xgb.DMatrix xgb.create.features xgb.importance xgb.dump
#' @importFrom Matrix t colSums
#' @importFrom methods as
#' @importFrom stats na.omit predict
#' @importFrom utils head tail
#' @export
dropSplit <- function(counts, score_cutoff = 0.9, GiniThreshold = NULL,
                      Uncertain_downsample = FALSE, Uncertain_downsample_times = 6, predict_Uncertain_only = TRUE,
                      Cell_rank = NULL, Uncertain_rank = NULL, Empty_rank = NULL,
                      modelOpt = FALSE, xgb_params = NULL, xgb_nrounds = 20, xgb_early_stopping_rounds = 3, xgb_thread = 8,
                      bounds = list(),
                      xgb_nfold = 5, xgb_metric = "auc",
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
    warning("'counts' has cells that nCount <=0. These cells will be removed in the following steps.",
      immediate. = TRUE
    )
    meta_info <- meta_info[meta_info$nCount > 0, ]
  }
  raw_cell_order <- rownames(meta_info)
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  counts <- counts[, rownames(meta_info)]

  nCount_inflection <- tail(find_inflection(meta_info$nCount)$index, 1)
  nFeature_inflection <- tail(find_inflection(meta_info$nFeature)$index, 1)
  inflection <- max(nCount_inflection, nFeature_inflection)
  meta_info$RankSE <- (meta_info$nCount_rank - meta_info$nFeature_rank)^2
  x <- meta_info$RankSE
  for (t in seq_len(3)) {
    x <- runMean(x, n = 100)
    x[is.na(x)] <- na.omit(x)[1]
  }
  meta_info$RankMSE <- x
  if (is.null(Cell_rank)) {
    inflection_left <- inflection - inflection * 0.4
    inflection_right <- inflection + inflection * 0.2
    Cell_rank <- inflection_left + which.min(meta_info$RankMSE[(inflection_left + 1):inflection_right]) - 1
  }
  if (is.null(Uncertain_rank)) {
    Uncertain_rank <- Cell_rank + which.max(meta_info$RankMSE[Cell_rank:(3 * Cell_rank)]) - 1
  }
  if (is.null(Empty_rank)) {
    Empty_rank <- Uncertain_rank + which.min(meta_info$RankMSE[Uncertain_rank:(5 * Uncertain_rank)]) - 1
  }

  Uncertain_count <- meta_info$nCount[Uncertain_rank]
  Cell_count <- meta_info$nCount[Cell_rank]
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
  # ggplot(meta_info, aes(x = log10(1:nrow(meta_info)), y = log10(RankMSE))) +
  #   geom_point(alpha = 0.1) +
  #   geom_vline(xintercept = log10(c(Cell_rank, Uncertain_rank, Empty_rank)))

  if (ncol(Empty_counts) > 100000) {
    warning("Too many Empty droplets. Only take the top 100000 droplets in the following steps.",
      immediate. = TRUE
    )
    Empty_counts <- Empty_counts[, 1:100000]
  }
  Uncertain_nCount <- Matrix::colSums(Uncertain_counts)
  Empty_nCount <- Matrix::colSums(Empty_counts)

  message(">>> Downsample pre-defined Cell droplets to a depth similar to the 'Empty': Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount))
  i <- sample(x = 1:ncol(Cell_counts), size = ncol(Empty_counts), replace = TRUE)
  Sim_Cell_counts <- Cell_counts[, i]
  colnames(Sim_Cell_counts) <- paste0("Sim-", 1:ncol(Sim_Cell_counts))
  Sim_Cell_nCount <- Matrix::colSums(Sim_Cell_counts)
  Sim_Cell_nCount_assign <- sample(unique(Empty_nCount), ncol(Sim_Cell_counts), replace = TRUE)
  Sim_Cell_counts <- downsampleMatrix(x = Sim_Cell_counts, prop = Sim_Cell_nCount_assign / Sim_Cell_nCount, bycol = TRUE)

  if (isTRUE(Uncertain_downsample)) {
    message(">>> Downsample pre-defined Uncertain droplets to a depth around 'Empty' ...")
    Sim_Uncertain_counts_list <- list()
    for (i in seq_len(Uncertain_downsample_times)) {
      x <- Uncertain_counts
      colnames(x) <- paste0("Sim", i, "-", colnames(x))
      Sim_Uncertain_nCount <- Uncertain_nCount
      Sim_Uncertain_nCount_assign <- sample(unique(Empty_nCount), ncol(x), replace = TRUE)
      x <- downsampleMatrix(x = x, prop = Sim_Uncertain_nCount_assign / Sim_Uncertain_nCount, bycol = TRUE)
      Sim_Uncertain_counts_list[[paste0("Sim-", i)]] <- x
    }
    Sim_Uncertain_counts <- Reduce(function(x, y) cbind(x, y), Sim_Uncertain_counts_list)
    comb_counts <- cbind(Cell_counts, Sim_Cell_counts, Sim_Uncertain_counts, Uncertain_counts, Empty_counts)
    group <- c(
      rep("Cell", ncol(Cell_counts)),
      rep("Sim_Cell", ncol(Sim_Cell_counts)),
      rep("Sim_Uncertain", ncol(Sim_Uncertain_counts)),
      rep("Uncertain", ncol(Uncertain_counts)),
      rep("Empty", ncol(Empty_counts))
    )
  } else {
    message(">>> Use raw Uncertain droplets...")
    Uncertain_downsample_times <- 1
    Sim_Uncertain_counts <- Uncertain_counts
    comb_counts <- cbind(Cell_counts, Sim_Cell_counts, Uncertain_counts, Empty_counts)
    group <- c(
      rep("Cell", ncol(Cell_counts)),
      rep("Sim_Cell", ncol(Sim_Cell_counts)),
      rep("Uncertain", ncol(Uncertain_counts)),
      rep("Empty", ncol(Empty_counts))
    )
  }

  message(">>> Calculate other cell metrics for the droplets to be trained...")
  comb_CellEntropy <- CellEntropy(comb_counts)
  comb_EntropyRate <- comb_CellEntropy / maxCellEntropy(comb_counts)
  comb_EntropyRate[is.na(comb_EntropyRate)] <- 1
  comb_Gini <- CellGini(comb_counts, normalize = T)
  if (is.null(GiniThreshold)) {
    GiniThreshold <- min(quantile(comb_Gini[colnames(Cell_counts)], 0.01), 0.99)
  }
  comb_GiniScore <- GiniScore(
    x = comb_Gini, GiniThreshold = GiniThreshold,
    group = group
  )
  MTgene <- grep(x = rownames(comb_counts), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T, value = TRUE)
  RPgene <- grep(x = rownames(comb_counts), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T, value = TRUE)
  comb_MTprop <- Matrix::colSums(comb_counts[MTgene, ]) / Matrix::colSums(comb_counts)
  comb_RPprop <- Matrix::colSums(comb_counts[RPgene, ]) / Matrix::colSums(comb_counts)
  dat_other <- cbind(
    CellEntropy = comb_CellEntropy,
    CellEntropyRate = comb_EntropyRate,
    MTprop = comb_MTprop,
    RPprop = comb_RPprop,
    CellGini = comb_Gini,
    GiniScore = comb_GiniScore
  )
  rownames(dat_other) <- colnames(comb_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_other)
  ] <- dat_other[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]

  norm_counts <- Matrix::t(Matrix::t(comb_counts) / Matrix::colSums(comb_counts))
  dat <- cbind(Matrix::t(norm_counts), dat_other[, !colnames(dat_other) %in% c("CellEntropy", "CellGini", "GiniScore")])
  ini_train <- dat[c(colnames(Cell_counts), colnames(Sim_Cell_counts), colnames(Empty_counts)), ]
  ini_train_label <- c(rep(1, ncol(Cell_counts) + ncol(Sim_Cell_counts)), rep(0, ncol(Empty_counts)))
  to_predict <- dat[c(colnames(Cell_counts), colnames(Sim_Uncertain_counts), colnames(Empty_counts)), ]

  if (isTRUE(modelOpt)) {
    opt <- xgbOptimization(
      dat = ini_train, dat_label = ini_train_label,
      bounds = bounds,
      xgb_nfold = xgb_nfold, xgb_nrounds = xgb_nrounds, xgb_early_stopping_rounds = xgb_early_stopping_rounds, xgb_metric = xgb_metric, xgb_thread = xgb_thread,
      opt_initPoints = opt_initPoints, opt_itersn = opt_itersn, opt_thread = opt_thread, ...
    )
    xgb_params <- c(opt$BestPars,
      eval_metric = "logloss",
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
      eval_metric = "logloss",
      eval_metric = "error",
      eval_metric = "auc",
      objective = "binary:logistic",
      nthread = xgb_thread
    )
  }
  message(">>> Prepare the XGBoost model with pre-defined classification...")
  # for (i in seq_len(10)) {
  # new_empty <- to_predict[rownames(meta_info)[meta_info$preDefinedClass=="Uncertain"&meta_info$dropSplitClass=="Empty"],]
  # message("Add ",nrow(new_empty)," new Empty droplets into training data from Uncertain droplets.")
  # train <- rbind(ini_train,new_empty)
  # train_label <- c(ini_train_label,rep(0,nrow(new_empty)))

  train <- ini_train
  train_label <- ini_train_label
  xgb <- xgboost(
    data = xgb.DMatrix(data = train, label = train_label),
    nrounds = xgb_nrounds,
    early_stopping_rounds = xgb_early_stopping_rounds,
    params = xgb_params
  )
  er <- tail(xgb$evaluation_log$train_error, 1)

  # new_train <- xgb.create.features(model = xgb, data=train)
  # new_to_predict <- xgb.create.features(model = xgb, data = to_predict)
  # new_xgb <- xgboost(
  #   data = xgb.DMatrix(data = new_train, label = train_label),
  #   nrounds = xgb_nrounds,
  #   early_stopping_rounds = xgb_early_stopping_rounds,
  #   params = xgb_params
  # )
  # if (nrow(to_predict) == 0) {
  #   message("No Uncertain droplets.")
  #   XGBoostScore <- numeric(0)
  # } else {
  #   XGBoostScore <- predict(xgb, to_predict)
  # }
  # XGBoostScore <- c(rep(1, ncol(Cell_counts)), XGBoostScore)

  XGBoostScore <- predict(xgb, to_predict)
  XGBoostScore <- c(
    XGBoostScore[1:ncol(Cell_counts)],
    rowMeans(matrix(XGBoostScore[(ncol(Cell_counts) + 1):(ncol(Cell_counts) + ncol(Sim_Uncertain_counts))], ncol = Uncertain_downsample_times)),
    XGBoostScore[(ncol(Cell_counts) + ncol(Sim_Uncertain_counts) + 1):length(XGBoostScore)]
  )
  if (isTRUE(predict_Uncertain_only)) {
    XGBoostScore[1:ncol(Cell_counts)] <- 1
    multiscore <- (XGBoostScore + ifelse(XGBoostScore > 0.5, 1, -1) * meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "GiniScore"] + ifelse(XGBoostScore > 0.5, 0, 1)) / 2
    score <- ifelse(XGBoostScore > 0.5, pmin(multiscore, XGBoostScore), pmax(multiscore, XGBoostScore))
    score[(ncol(Cell_counts) + ncol(Uncertain_counts) + 1):length(score)][score[(ncol(Cell_counts) + ncol(Uncertain_counts) + 1):length(score)] > score_cutoff] <- 0.5
  } else {
    multiscore <- (XGBoostScore + ifelse(XGBoostScore > 0.5, 1, -1) * meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "GiniScore"] + ifelse(XGBoostScore > 0.5, 0, 1)) / 2
    score <- ifelse(XGBoostScore > 0.5, pmin(multiscore, XGBoostScore), pmax(multiscore, XGBoostScore))
  }

  meta_info[, "preDefinedClass"] <- "Discarded"
  meta_info[colnames(Cell_counts), "preDefinedClass"] <- "Cell"
  meta_info[colnames(Uncertain_counts), "preDefinedClass"] <- "Uncertain"
  meta_info[colnames(Empty_counts), "preDefinedClass"] <- "Empty"
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "XGBoostScore"] <- XGBoostScore
  meta_info[, "dropSplitClass"] <- meta_info[, "preDefinedClass"]
  meta_info[, "dropSplitScore"] <- 0
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "dropSplitScore"] <- score
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "dropSplitClass"] <- ifelse(
    score > score_cutoff, "Cell", ifelse(score < 1 - score_cutoff, "Empty", "Uncertain")
  )

  rescure <- which(meta_info[, "preDefinedClass"] %in% c("Uncertain", "Empty") & meta_info[, "dropSplitClass"] == "Cell")
  rescure_score <- meta_info[rescure, "dropSplitScore"]
  er_rate <- min((1 - score_cutoff) / (1 - er * 2), 1)
  drop <- ceiling(length(rescure) * er_rate)
  drop_index <- rescure[order(rescure_score)[1:drop]]
  message(
    "\n>>> Control the rate of false positives:",
    "\n... Number of new defined Cell from Uncertain or Empty:", length(rescure),
    "\n... Estimated error rate:", round(er_rate, digits = 6),
    "\n... Estimated error number:", drop,
    "\n... Error droplets mean dropSplitScore:", round(mean(meta_info[drop_index, "dropSplitScore"]), digits = 6)
  )
  meta_info[drop_index, "dropSplitScore"] <- 0.5
  meta_info[drop_index, "dropSplitClass"] <- "Uncertain"

  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
  meta_info <- meta_info[raw_cell_order, ]

  message(
    "\n>>> Summary of dropSplit-defined droplet:",
    "\n... Number of Cell: ", sum(meta_info$dropSplitClass == "Cell"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Cell", "nCount"]),
    "\n... Number of Uncertain: ", sum(meta_info$dropSplitClass == "Uncertain"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Uncertain", "nCount"]),
    "\n... Number of Empty: ", sum(meta_info$dropSplitClass == "Empty"), "  Minimum nCounts: ", min(meta_info[meta_info$dropSplitClass == "Empty", "nCount"]),
    "\n... Number of Discarded: ", sum(meta_info$dropSplitClass == "Discarded"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Discarded") > 1, min(meta_info[meta_info$dropSplitClass == "Discarded", "nCount"]), min(meta_info[, "nCount"])),
    "\n>>> Pre-defined as 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell"),
    "\n... Mean GiniScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "GiniScore"]), 3),
    "\n... Mean XGBoostScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "XGBoostScore"]), 3),
    "\n... Mean dropSplitScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "dropSplitScore"]), 3),
    "\n>>> Pre-defined as 'Uncertain' switch to 'Cell': ", sum(meta_info$preDefinedClass == "Uncertain" & meta_info$dropSplitClass == "Cell"),
    "\n... Mean GiniScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass == "Cell", "GiniScore"]), 3),
    "\n... Mean XGBoostScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
    "\n... Mean dropSplitScore:", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3)
  )
  # }

  importance_matrix <- as.data.frame(xgb.importance(model = xgb))
  tree <- xgb.dump(xgb, with_stats = TRUE)
  result <- list(
    meta_info = meta_info,
    train = train,
    train_label = train_label,
    to_predict = to_predict,
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
