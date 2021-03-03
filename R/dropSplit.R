#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of three parts: quality control, model construction and droplet classification, summarizing features.
#' dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent samples.
#' @param score_cutoff A cutoff value of \code{dropSplitScore} to determine if a droplet is cell-containing or empty. Default is 0.9, i.e. a droplet with \code{dropSplitScore}>0.9 will be classified as \code{Cell} and a with \code{dropSplitScore}<0.1  will be classified as \code{Empty}.
#' @param Gini_control Whether to control cell quality by CellGini. Default is \code{TRUE}.
#' @param Gini_threshold A value used in \code{\link{Score}} function for CellGini metric. The higher, the more conservative and will get a lower number of cells. Default is automatic.
#' @param Uncertain_downsample Whether to use a downsampled Uncertain droplets for predicting. Default is FALSE.
#' @param Uncertain_downsample_times Number of downsample times for each Uncertain droplet. \code{dropSplitScore} of downsampled droplets from the same Uncertain droplet will be averaged. Default is 6.
#' @param predict_Uncertain_only Whether to predict only the Uncertain droplets. Default is \code{TRUE}.
#' @param remove_FP_by metric used to remove the estimated false positives by. Must be one of \code{nCount}, \code{nFeature}, \code{CellEntropy}, \code{CellEfficiency}, \code{dropSplitScore}. Default is \code{dropSplitScore}.
#' @param Cell_rank,Uncertain_rank,Empty_rank Custom Rank value to mark the droplets as Cell, Uncertain and Empty labels for the data to be trained. Default is automatic. But useful when the default value is considered to be wrong from the RankMSE plot.
#' @param modelOpt Whether to optimize the model using \code{\link{xgbOptimization}}. Will take long time for large datasets. If \code{TRUE}, will overwrite the parameters list in \code{xgb_params}.
#' @param xgb_params The \code{list} of XGBoost parameters.
#' @inheritParams xgbOptimization
#' @param ... Other arguments passed to \code{\link{xgbOptimization}}.
#'
#' @return A list of seven objects:
#' \describe{
#' \item{meta_info}{A \code{DataFrame} object of evaluation metrics to be used in dropSplit and the final droplet classification.}
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
#' @importFrom methods hasArg
#' @importFrom S4Vectors DataFrame
#' @export
dropSplit <- function(counts, score_cutoff = 0.9, Gini_control = TRUE, Gini_threshold = 0.99,
                      Uncertain_downsample = FALSE, Uncertain_downsample_times = 6, predict_Uncertain_only = TRUE, remove_FP_by = "dropSplitScore",
                      Cell_rank = NULL, Uncertain_rank = NULL, Empty_rank = NULL,
                      modelOpt = FALSE, xgb_params = NULL, xgb_nrounds = 20, xgb_early_stopping_rounds = 3, xgb_thread = 8,
                      bounds = list(),
                      xgb_nfold = 5, xgb_metric = "auc",
                      opt_initPoints = length(bounds) + 1, opt_itersn = 10, opt_thread = 1, ...) {
  if (!hasArg(counts)) {
    stop("Parameter 'counts' not found.")
  }
  if (!class(counts) %in% c("matrix", "dgCMatrix", "dgTMatrix")) {
    stop("'counts' must be a dense matrix or a sparse matrix object.")
  }
  if (length(colnames(counts)) != ncol(counts) | length(rownames(counts)) != nrow(counts)) {
    stop("'counts' matrix must have both row(feature) names and column(cell) names.")
  }
  if (!is.null(score_cutoff)) {
    if (score_cutoff < 0 | score_cutoff > 1) {
      stop("'score_cutoff' must be between 0 and 1.")
    }
  }
  if (!is.null(Gini_threshold)) {
    if (Gini_threshold < 0 | Gini_threshold > 1) {
      stop("'Gini_threshold' must be between 0 and 1.")
    }
  }
  if (!remove_FP_by %in% c("nCount", "nFeature", "CellEntropy", "CellEfficiency", "CellRedundancy", "dropSplitScore") | length(remove_FP_by) > 1) {
    stop("'remove_FP_by' must be one of nCount, nFeature, CellEntropy, CellEfficiency, CellRedundancy, dropSplitScore.")
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
    warning("'counts' has droplets that nCount <=0. These droplets are removed in the following steps.",
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
    inflection_left <- inflection - inflection * 0.3
    inflection_right <- inflection + inflection * 0.2
    Cell_rank <- inflection_left + which.min(meta_info$RankMSE[(inflection_left + 1):inflection_right]) - 1
  }
  if (is.null(Uncertain_rank)) {
    Uncertain_rank <- Cell_rank + which.max(meta_info$RankMSE[Cell_rank:min(3 * Cell_rank, nrow(meta_info))]) - 1
  }
  if (is.null(Empty_rank)) {
    Empty_rank <- Uncertain_rank + which.min(meta_info$RankMSE[Uncertain_rank:min(10 * Cell_rank, nrow(meta_info))]) - 1
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

  message(">>> Calculate various cell metrics for the droplets...")
  final_counts <- cbind(Cell_counts, Uncertain_counts, Empty_counts)
  final_Gini <- CellGini(final_counts, normalize = TRUE)
  if (is.null(Gini_threshold)) {
    Gini_threshold <- min(0.99, quantile(final_Gini[colnames(Cell_counts)], 0.01))
  }
  final_CellGiniScore <- Score(
    x = final_Gini, threshold = Gini_threshold,
    group = c(
      rep("1", ncol(Cell_counts)),
      rep("2", ncol(Uncertain_counts)),
      rep("3", ncol(Empty_counts))
    )
  )
  final_CellEntropy <- CellEntropy(final_counts)
  final_maxCellEntropy <- maxCellEntropy(final_counts)
  final_CellEfficiency <- final_CellEntropy / final_maxCellEntropy
  final_CellEfficiency[is.na(final_CellEfficiency)] <- 1
  final_CellRedundancy <- 1 - final_CellEfficiency
  final_CellRedundancy[final_CellRedundancy < 0] <- 0
  dat_metrics <- cbind(
    CellEntropy = final_CellEntropy,
    CellEfficiency = final_CellEfficiency,
    CellRedundancy = final_CellRedundancy,
    CellGini = final_Gini,
    CellGiniScore = final_CellGiniScore
  )
  rownames(dat_metrics) <- colnames(final_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_metrics)
  ] <- dat_metrics[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]

  cell_Gini <- final_Gini[colnames(Cell_counts)]
  cell_drop <- names(cell_Gini)[cell_Gini < Gini_threshold]
  cell_use <- names(cell_Gini)[cell_Gini >= Gini_threshold]
  if (length(cell_drop) > 0) {
    message(
      "*** There are ", length(cell_drop), " pre-defined droplets with lower CellGini value than Gini_threshold:", round(Gini_threshold, 3), "\n",
      "*** These ", length(cell_drop), " droplets are converted to 'Uncertain'."
    )
    Drop_counts <- Cell_counts[, cell_drop]
    Cell_counts <- Cell_counts[, cell_use]
    Uncertain_counts <- cbind(Drop_counts, Uncertain_counts)
  }

  message(
    ">>> Downsample pre-defined Cell droplets to a depth similar to the 'Empty'\n",
    "... nCount in 'Empty' droplets: Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount)
  )
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
  } else {
    message(">>> Use raw Uncertain droplets...")
    Uncertain_downsample_times <- 1
    Sim_Uncertain_counts <- Uncertain_counts
    comb_counts <- cbind(Cell_counts, Sim_Cell_counts, Uncertain_counts, Empty_counts)
  }

  message(">>> Calculate MTprop and RPprop for the droplets to be trained...")
  MTgene <- grep(x = rownames(comb_counts), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T, value = TRUE)
  RPgene <- grep(x = rownames(comb_counts), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T, value = TRUE)
  comb_MTprop <- Matrix::colSums(comb_counts[MTgene, ]) / Matrix::colSums(comb_counts)
  comb_RPprop <- Matrix::colSums(comb_counts[RPgene, ]) / Matrix::colSums(comb_counts)
  dat_Features <- cbind(
    MTprop = comb_MTprop,
    RPprop = comb_RPprop
  )
  rownames(dat_Features) <- colnames(comb_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_Features)
  ] <- dat_Features[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]


  message(">>> Merge new features into train data...")
  norm_counts <- Matrix::t(Matrix::t(comb_counts) / Matrix::colSums(comb_counts))
  dat <- cbind(Matrix::t(norm_counts), dat_Features[, colnames(dat_Features)])
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
  message(">>> Construct the XGBoost model with pre-defined classification...")
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
  names(XGBoostScore) <- c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts))

  if (isTRUE(predict_Uncertain_only)) {
    XGBoostScore[colnames(Cell_counts)] <- 1
    XGBoostScore[colnames(Empty_counts)][XGBoostScore[colnames(Empty_counts)] > 0.5] <- 0.5
  }

  multiscore <- XGBoostScore
  if (isTRUE(Gini_control)) {
    raw_score <- multiscore
    multiscore <- (multiscore + ifelse(multiscore > 0.5, 1, -1) * meta_info[names(multiscore), "CellGiniScore"] + ifelse(multiscore > 0.5, 0, 1)) / 2
    multiscore <- ifelse(raw_score > 0.5, pmin(raw_score, multiscore), pmax(raw_score, multiscore))
  }

  meta_info[, "preDefinedClass"] <- "Discarded"
  meta_info[colnames(Cell_counts), "preDefinedClass"] <- "Cell"
  meta_info[colnames(Uncertain_counts), "preDefinedClass"] <- "Uncertain"
  meta_info[colnames(Empty_counts), "preDefinedClass"] <- "Empty"
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "XGBoostScore"] <- XGBoostScore
  meta_info[, "dropSplitClass"] <- meta_info[, "preDefinedClass"]
  meta_info[, "dropSplitScore"] <- 0
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "dropSplitScore"] <- multiscore
  meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "dropSplitClass"] <- ifelse(
    multiscore > score_cutoff, "Cell", ifelse(multiscore < 1 - score_cutoff, "Empty", "Uncertain")
  )

  if (remove_FP_by %in% c("nCount", "nFeature", "CellEntropy", "CellRedundancy", "dropSplitScore")) {
    decreasing <- FALSE
  }
  if (remove_FP_by %in% c("CellEfficiency")) {
    decreasing <- TRUE
  }
  rescure <- which(meta_info[, "preDefinedClass"] %in% c("Uncertain", "Empty") & meta_info[, "dropSplitClass"] == "Cell")
  rescure_score <- meta_info[rescure, remove_FP_by]
  er_rate <- min((1 - score_cutoff) / (1 - er * 2), 1)
  drop <- ceiling(length(rescure) * er_rate)
  if (drop > 0) {
    drop_index <- rescure[order(rescure_score, decreasing = decreasing)[1:drop]]
    drop_score <- round(mean(meta_info[drop_index, remove_FP_by]), digits = 3)
    meta_info[drop_index, "dropSplitScore"] <- 0.5
    meta_info[drop_index, "dropSplitClass"] <- "Uncertain"
  } else {
    drop_score <- NA
  }
  message(
    "\n>>> Control the rate of false positives",
    "\n... Number of new defined Cell from Uncertain or Empty: ", length(rescure),
    "\n... Estimated error rate: ", round(er_rate, digits = 3),
    "\n... Estimated error number: ", drop,
    "\n... Estimated error droplets mean ", remove_FP_by, ": ", drop_score,
    "\n*** The dropSplitScore and dropSplitClass for these ", drop, " droplets are converted to 0.5 and 'Uncertain'"
  )

  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
  meta_info <- meta_info[raw_cell_order, ]

  message(
    "\n>>> Summary of dropSplit-defined droplet",
    "\n... Number of 'Cell': ", sum(meta_info$dropSplitClass == "Cell"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Cell") > 1, min(meta_info[meta_info$dropSplitClass == "Cell", "nCount"]), 1),
    "\n... Number of 'Uncertain': ", sum(meta_info$dropSplitClass == "Uncertain"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Uncertain") > 1, min(meta_info[meta_info$dropSplitClass == "Uncertain", "nCount"]), 1),
    "\n... Number of 'Empty': ", sum(meta_info$dropSplitClass == "Empty"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Empty") > 1, min(meta_info[meta_info$dropSplitClass == "Empty", "nCount"]), 1),
    "\n... Number of 'Discarded': ", sum(meta_info$dropSplitClass == "Discarded"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Discarded") > 1, min(meta_info[meta_info$dropSplitClass == "Discarded", "nCount"]), 1),
    "\n>>> Pre-defined as 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell"),
    "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "XGBoostScore"]), 3),
    "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "CellGiniScore"]), 3),
    "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$preDefinedClass == "Cell" & meta_info$dropSplitClass != "Cell", "dropSplitScore"]), 3),
    "\n>>> Pre-defined as 'Uncertain' or 'Empty' switch to 'Cell': ", sum(meta_info$preDefinedClass != "Cell" & meta_info$dropSplitClass == "Cell"),
    "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$preDefinedClass != "Cell" & meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
    "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$preDefinedClass != "Cell" & meta_info$dropSplitClass == "Cell", "CellGiniScore"]), 3),
    "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$preDefinedClass != "Cell" & meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3)
  )
  # }

  importance_matrix <- as.data.frame(xgb.importance(model = xgb))
  tree <- xgb.dump(xgb, with_stats = TRUE)
  result <- list(
    meta_info = DataFrame(meta_info),
    train = train,
    train_label = train_label,
    to_predict = to_predict,
    model = xgb,
    importance_matrix = importance_matrix,
    tree = tree
  )

  return(result)
}
