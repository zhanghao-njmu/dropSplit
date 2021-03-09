#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of three parts: quality control, model construction and droplet classification, summarizing features.
#' dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent droplets and rows represent features.
#' @param do_plot Whether to plot during the cellcalling. Default is \code{TRUE}.
#' @param cell_score A cutoff value of \code{dropSplitScore} to determine if a droplet is cell-containing. Range between 0.5 and 1. Default is 0.8.
#' @param empty_score A cutoff value of \code{dropSplitScore} to determine if a droplet is empty. Range between 0 and 0.5. Note, a reasonable will also improve the accuracy for 'Cell' identification and this can be checked by \code{\link{QCplot}}. Default is 0.2.
#' @param CE_ratio Ratio value between down-sampled 'Cells' and 'Empty' droplets. The actual value will be slightly higher than the set. \code{CE_ratio>1} can help with 'Cell' identification. Default is 2.
#' @param max_iter An integer specifying the number of iterations to use to rebuild the model with new defined droplets. Default is 10.
#' @param min_error The minimum train error value to be achieved by the model. Default is 0.002.
#' @param min_improve Minimal improvement of the model. Default is 0.002.
#' @inheritParams RankMSE
#' @param Gini_control Whether to control cell quality by CellGini. Default is \code{TRUE}.
#' @param Gini_threshold A value used in \code{\link{Score}} function for CellGini metric. The higher, the more conservative and will get a lower number of cells. Default is automatic.
#' @param Cell_rank,Uncertain_rank,Empty_rank Custom Rank value to mark the droplets as Cell, Uncertain and Empty labels for the data to be trained. Default is automatic. But useful when the default value is considered to be wrong from the RankMSE plot.
#' @param Uncertain_downsample_times Number of downsample times for each Uncertain droplet. \code{XGBoostScore} of down-sampled droplets from the same Uncertain droplet will be averaged and perform \code{\link[stats]{t.test}} with \code{cell_score}. Default is 10.
#' @param cell_score_FDR FDR cutoff for droplets that predicted as 'Cell' from pre-defined 'Uncertain' or pre-defined 'Cell' droplets. Note, statistic tests and the FDR control only performed on the difference between averaged \code{XGBoostScore} and cell_score. It is not a false 'Cell' rate. Default is 0.05.
#' @param preCell_mask logical; Whether to mask pre-defined 'Cell' droplets when prediction. If \code{TRUE}, XGBoostScore for all droplets pre-defined as 'Cell' will be set to 1; Default is automatic.
#' @param preEmpty_mask logical; Whether to mask pre-defined 'Empty' droplets when prediction. There is a little different with parameter \code{preCell_mask}. If \code{TRUE}, XGBoostScore will not change, but the final classification will not be 'Cell' in any case. Default is \code{TRUE}.
#' @param xgb_params The \code{list} of XGBoost parameters.
#' @param modelOpt Whether to optimize the model using \code{\link{xgbOptimization}}. Will take long time for large datasets. If \code{TRUE}, will overwrite the parameters list in \code{xgb_params}. The following parameters are only used in \code{\link{xgbOptimization}}.
#' @inheritParams xgbOptimization
#' @param seed Random seed used in simulation. Default is 0.
#' @param ... Other arguments passed to \code{\link{xgbOptimization}}.
#'
#' @return A list of six objects:
#' \describe{
#' \item{meta_info}{A \code{DataFrame} object of evaluation metrics to be used in dropSplit and classification task. Important columns in the \code{meta_info}:
#' \itemize{
#' \item \code{preDefinedClass}
#' \item \code{CellGini}
#' \item \code{CellGiniScore}
#' \item \code{XGBoostScore}
#' \item \code{pvalue}
#' \item \code{FDR}
#' \item \code{dropSplitClass}
#' \item \code{dropSplitScore}
#' }}
#' \item{train}{The dataset trained in the final XGBoost model. It consists of two pre-defined droplets: Cell(raw + simulated) and Empty.}
#' \item{train_label}{Labels for the \code{train}. 0 represents 'Empty', 1 represents 'Cell'.}
#' \item{to_predict}{The dataset that to be predicted. It consists of all three pre-defined droplets: Cell(raw + simulated), Uncertain(simulated) and Empty.}
#' \item{model}{The XGBoost model used in dropSplit for classification.}
#' \item{importance_matrix}{A \code{data.frame} of feature importances in the classification model.}
#' }
#'
#' @examples
#' library(dropSplit)
#' # Simulate a counts matrix including 20000 empty droplets, 2000 large cells and 200 small cells.
#' simple_counts <- simSimpleCounts()
#' true <- strsplit(colnames(simple_counts), "-")
#' true <- as.data.frame(Reduce(function(x, y) rbind(x, y), true))
#' colnames(true) <- c("label", "Type", "Cluster", "Cell")
#' rownames(true) <- colnames(simple_counts)
#' true_label <- true$label
#' table(true_label)
#'
#' ## DropSplit ---------------------------------------------------------------
#' result <- dropSplit(simple_counts)
#' qc <- QCplot(result$meta_info)
#' dropSplitClass <- result$meta_info$dropSplitClass
#' table(true_label, result$meta_info$dropSplitClass)
#'
#' # compare with the true labels
#' result$meta_info$true_label <- true_label
#' qc_true <- QCplot(result$meta_info, colorBy = "true_label")
#' qc_true$CellEntropy$Merge
#'
#' # QC plot using all metrics
#' qc <- QCplot(result$meta_info)
#' qc$RankMSE$Merge
#' qc$CellEntropy$Merge
#' qc$CellEfficiency$Merge
#'
#' # Feature importance plot
#' fp <- ImportancePlot(result$meta_info, result$train, result$importance_matrix, top_n = 20)
#' fp$Importance
#' fp$preDefinedClassExp
#' fp$dropSplitClassExp
#' @importFrom TTR runMean
#' @importFrom DropletUtils downsampleMatrix
#' @importFrom xgboost xgboost xgb.DMatrix xgb.create.features xgb.importance xgb.dump
#' @importFrom Matrix t colSums
#' @importFrom methods as
#' @importFrom stats na.omit predict median quantile t.test p.adjust
#' @importFrom utils head tail
#' @importFrom methods hasArg
#' @importFrom S4Vectors DataFrame
#' @importFrom ggplot2 ggplot aes geom_point geom_vline
#' @export
dropSplit <- function(counts, do_plot = TRUE, cell_score = 0.8, empty_score = 0.2, CE_ratio = 2,
                      smooth_num = 5, smooth_window = 100, tolerance = 0.5,
                      Gini_control = TRUE, Gini_threshold = NULL,
                      Cell_rank = NULL, Uncertain_rank = NULL, Empty_rank = NULL,
                      Uncertain_downsample_times = 10, cell_score_FDR = 0.05,
                      max_iter = 10, min_error = 0.002, min_improve = 0.002,
                      preCell_mask = NULL, preEmpty_mask = TRUE,
                      xgb_params = NULL, xgb_nrounds = 20, xgb_thread = 8,
                      modelOpt = FALSE, seed = 0, ...) {
  start <- Sys.time()
  if (!hasArg(counts)) {
    stop("Parameter 'counts' not found.")
  }
  if (!class(counts) %in% c("matrix", "dgCMatrix", "dgTMatrix")) {
    stop("'counts' must be a dense matrix or a sparse matrix object.")
  }
  if (length(colnames(counts)) != ncol(counts) | length(rownames(counts)) != nrow(counts)) {
    stop("'counts' matrix must have both row(feature) names and column(cell) names.")
  }
  if (!is.null(cell_score)) {
    if (cell_score < 0.5 | cell_score > 1) {
      stop("'cell_score' must be between 0.5 and 1.")
    }
  }
  if (!is.null(Gini_threshold)) {
    if (Gini_threshold < 0 | Gini_threshold > 1) {
      stop("'Gini_threshold' must be between 0 and 1.")
    }
  }
  if (Uncertain_downsample_times < 3) {
    stop("'Uncertain_downsample_times' must be larger than 3.")
  }
  if (class(counts) == "matrix") {
    counts <- as(counts, "dgCMatrix")
  }

  set.seed(seed)

  message(">>> Start to define the credible cell-containing droplet...")
  meta_info <- data.frame(row.names = colnames(counts))
  meta_info$nCount <- Matrix::colSums(counts)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  meta_info$nFeature <- Matrix::colSums(counts > 0)
  meta_info$nFeature_rank <- rank(-(meta_info$nFeature))

  if (min(meta_info$nCount) <= 0) {
    stop("'counts' has droplets that nCount<=0.")
  }
  raw_cell_order <- rownames(meta_info)
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  counts <- counts[, rownames(meta_info)]

  out <- RankMSE(meta_info = meta_info, smooth_num = smooth_num, smooth_window = smooth_window, find_cell_count = TRUE, tolerance = tolerance)
  meta_info <- out$meta_info
  cell_count <- out$cell_count
  # qplot(meta_info$nCount_rank,log(meta_info$RankMSE))+xlim(1,20000)+geom_vline(xintercept = (Cell_rank))

  if (is.null(Cell_rank)) {
    Cell_rank <- max(meta_info$nCount_rank[meta_info$nCount > cell_count])
  }
  if (is.null(Uncertain_rank)) {
    Uncertain_rank <- Cell_rank + which.max(meta_info$RankMSE[Cell_rank:min(Cell_rank * 2, nrow(meta_info))])[1] - 1
  }
  if (is.null(Empty_rank)) {
    Empty_rank <- Uncertain_rank * 2 + which.min(meta_info$RankMSE[(Uncertain_rank * 2):min(Uncertain_rank * 10, nrow(meta_info))])[1] - 1
  }
  Cell_count <- meta_info$nCount[Cell_rank]
  Uncertain_count <- meta_info$nCount[Uncertain_rank]
  Empty_count <- max(meta_info$nCount[Empty_rank], 10)

  Cell_counts <- counts[, meta_info$nCount >= Cell_count]
  Uncertain_counts <- counts[, meta_info$nCount < Cell_count & meta_info$nCount >= Uncertain_count]
  Empty_counts <- counts[, meta_info$nCount < Uncertain_count & meta_info$nCount >= Empty_count]

  if (ncol(Uncertain_counts) == 0) {
    stop("No 'Uncertain' droplets detected. Please check the RankMSE curve with the pre-defined droplet cutoff. You may set custom cutoff values in the parameters manually.")
  }
  if (ncol(Empty_counts) == 0) {
    stop("No 'Empty' droplets detected. Please check the RankMSE curve with the pre-defined droplet cutoff. You may set custom cutoff values in the parameters manually.")
  }
  if (ncol(Empty_counts) < ncol(Uncertain_counts) | ncol(Empty_counts) < ncol(Cell_counts)) {
    stop("'Empty' droplets is fewer than 'Cell' or 'Uncertain'. Please check the RankMSE curve with the pre-defined droplet cutoff. You may set custom cutoff values in the parameters manually.")
  }

  message(">>> Calculate various cell metrics for the droplets...")
  final_counts <- cbind(Cell_counts, Uncertain_counts, Empty_counts)
  final_Gini <- CellGini(final_counts, normalize = TRUE)
  if (is.null(Gini_threshold)) {
    Gini_threshold <- (quantile(final_Gini[colnames(Cell_counts)], 0.1) + (1 - 2 * cell_score) * quantile(final_Gini[colnames(Cell_counts)], 0.9)) / (2 - 2 * cell_score)
    Gini_threshold <- ifelse(Gini_threshold < 0.95, 0.95, ifelse(Gini_threshold > 0.99, 0.99, Gini_threshold))
  }
  message("*** Gini_threshold was set to: ", round(Gini_threshold, 3))
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
      "*** There are ", length(cell_drop), " high RNA content droplets with low 'CellGini' values than the Gini_threshold:", round(Gini_threshold, 3), "\n",
      "*** These ", length(cell_drop), " droplets will pre-defined as 'Uncertain'."
    )
    Drop_counts <- Cell_counts[, cell_drop]
    if (is.numeric(Drop_counts)) {
      Drop_counts <- as.matrix(Drop_counts)
      colnames(Drop_counts) <- cell_drop
    }
    Cell_counts <- Cell_counts[, cell_use]
    Uncertain_counts <- cbind(as(Drop_counts, "dgCMatrix"), Uncertain_counts)
  }

  if (ncol(Empty_counts) > 50000) {
    warning("Too many 'Empty' droplets. Only take the top 50000 'Empty' droplets by nCount in the following steps.",
            immediate. = TRUE
    )
    Empty_counts <- Empty_counts[, 1:50000]
  }

  Cell_nCount <- Matrix::colSums(Cell_counts)
  Uncertain_nCount <- Matrix::colSums(Uncertain_counts)
  Empty_nCount <- Matrix::colSums(Empty_counts)

  meta_info[, "preDefinedClass"] <- "Discarded"
  meta_info[colnames(Cell_counts), "preDefinedClass"] <- "Cell"
  meta_info[colnames(Uncertain_counts), "preDefinedClass"] <- "Uncertain"
  meta_info[colnames(Empty_counts), "preDefinedClass"] <- "Empty"
  message(
    ">>> Summary of pre-defined droplets",
    "\n... Number of Cell: ", ncol(Cell_counts), "  Minimum nCounts: ", Cell_count,
    "\n... Number of Uncertain: ", ncol(Uncertain_counts), "  Minimum nCounts: ", Uncertain_count,
    "\n... Number of Empty: ", ncol(Empty_counts), "  Minimum nCounts: ", Empty_count,
    "\n... Number of Discarded: ", ncol(counts) - ncol(Cell_counts) - ncol(Uncertain_counts) - ncol(Empty_counts), "  Minimum nCounts: ", min(meta_info$nCount)
  )
  if (do_plot) {
    p <- RankMSEPlot(meta_info, colorBy = "preDefinedClass", splitBy = NULL, cell_stat_by = "preDefinedClass")
    suppressWarnings(print(p))
  }

  message(
    ">>> Downsample pre-defined 'Cell' droplets to a depth similar to the 'Empty' for training",
    "\n... nCount in 'Cell' droplets: Min=", min(Cell_nCount), " Median=", median(Cell_nCount), " Max=", max(Cell_nCount),
    "\n... nCount in 'Empty' droplets: Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount)
  )
  # a <- as.data.frame(Cell_nCount)
  # a$Cell <- rownames(a)
  # b <- as.data.frame(table(Cell_nCount))
  # c <- merge(x = a, y = b, by = "Cell_nCount")
  # rownames(c) <- c$Cell
  # c$Rank <- meta_info[c$Cell, "nCount_rank"]
  # c$prop <- c$Freq / sum(c$Freq)

  Cell_downsample_times <- ceiling(ncol(Empty_counts) / ncol(Cell_counts) * CE_ratio)
  Sim_Cell_counts <- Cell_counts[, rep(1:ncol(Cell_counts), Cell_downsample_times)]
  colnames(Sim_Cell_counts) <- paste0(rep(paste0("Sim", 1:Cell_downsample_times), each = ncol(Cell_counts)), "-", colnames(Sim_Cell_counts))
  Sim_Cell_nCount_assign <- sample(unique(Empty_nCount), ncol(Sim_Cell_counts), replace = TRUE)
  Sim_Cell_counts <- downsampleMatrix(x = Sim_Cell_counts, prop = Sim_Cell_nCount_assign / Matrix::colSums(Sim_Cell_counts), bycol = TRUE)
  Sim_Cell_nCount <- Matrix::colSums(Sim_Cell_counts)

  message(
    ">>> Downsample pre-defined 'Uncertain' droplets to a depth similar to the 'Empty' for prediction",
    "\n... nCount in 'Uncertain' droplets: Min=", min(Uncertain_nCount), " Median=", median(Uncertain_nCount), " Max=", max(Uncertain_nCount),
    "\n... nCount in 'Empty' droplets: Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount),
    "\n... Downsample times for each 'Uncertain' droplet: ", Uncertain_downsample_times
  )
  Sim_Uncertain_counts <- Uncertain_counts[, rep(1:ncol(Uncertain_counts), Uncertain_downsample_times)]
  colnames(Sim_Uncertain_counts) <- paste0(rep(paste0("Sim", 1:Uncertain_downsample_times), each = ncol(Uncertain_counts)), "-", colnames(Sim_Uncertain_counts))
  Sim_Uncertain_nCount_assign <- sample(unique(Empty_nCount), ncol(Sim_Uncertain_counts), replace = TRUE)
  Sim_Uncertain_counts <- downsampleMatrix(x = Sim_Uncertain_counts, prop = Sim_Uncertain_nCount_assign / Matrix::colSums(Sim_Uncertain_counts), bycol = TRUE)
  Sim_Uncertain_nCount <- Matrix::colSums(Sim_Uncertain_counts)

  message(">>> Calculate MTprop and RPprop for the droplets to be trained...")
  comb_counts <- cbind(Cell_counts, Sim_Cell_counts, Sim_Uncertain_counts, Uncertain_counts, Empty_counts)
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
  to_predict <- dat[c(colnames(Cell_counts), colnames(Sim_Cell_counts), colnames(Sim_Uncertain_counts), colnames(Empty_counts)), ]
  Sim_Uncertain_counts_rawname <- gsub(x = colnames(Sim_Uncertain_counts), pattern = "Sim\\d+-", replacement = "")

  if (isTRUE(modelOpt)) {
    opt <- xgbOptimization(
      dat = ini_train, dat_label = ini_train_label,
      xgb_nrounds = xgb_nrounds, xgb_thread = xgb_thread, ...
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
  message(">>> Construct the XGBoost model with pre-defined classification...")

  train_error <- 1
  k <- 0
  j <- 0
  if (max_iter == 1) {
    message("max_iter=1 but at least 2 iteration needed. max_iter is set to 2.")
    max_iter <- 2
  }
  while (k < max_iter) {
    k <- k + 1
    j <- j + 1
    message("\n  ============= Iteration: ", k, " =============  ")
    if (k == 1) {
      train <- ini_train
      train_label <- ini_train_label
    } else {
      empty_current <- rownames(meta_info)[meta_info$dropSplitClass == "Empty"]
      raw_empty <- colnames(Empty_counts)[colnames(Empty_counts) %in% empty_current]
      new_empty <- colnames(Sim_Uncertain_counts)[Sim_Uncertain_counts_rawname %in% empty_current]
      empty_update <- c(new_empty, raw_empty)
      if (length(empty_update) > ncol(Empty_counts)) {
        message("*** Number of 'Empty' droplets is too large. Take the top ", ncol(Empty_counts), " by nCount for training.")
        empty_update_nCount <- c(Sim_Uncertain_nCount[new_empty], Empty_nCount[raw_empty])
        empty_update_order <- order(empty_update_nCount, decreasing = TRUE)
        empty_update <- empty_update[empty_update_order[1:ncol(Empty_counts)]]
      }
      if (length(empty_update) < ncol(Cell_counts)) {
        message("*** Number of 'Empty' droplets is too small. Use the previous model for final classification.")
        break
      }
      message(
        "... Number of removed pre-defined 'Empty' droplets: ", ncol(Empty_counts) - length(raw_empty),
        "\n... Number of newly added 'Empty' droplets from 'Uncertain'(including the simulated): ", length(new_empty),
        "\n... Number of total 'Empty' droplets used for taining (including the simulated): ", length(empty_update)
      )
      train <- dat[c(colnames(Cell_counts), colnames(Sim_Cell_counts), empty_update), ]
      train_label <- c(rep(1, ncol(Cell_counts) + ncol(Sim_Cell_counts)), rep(0, length(empty_update)))
    }

    message(">>> Train data: 'Cell'=", sum(train_label == 1), "; 'Empty'=", sum(train_label == 0), "; 'C/E' Ratio=", round(sum(train_label == 1) / sum(train_label == 0), 3))
    xgb <- xgboost(
      data = xgb.DMatrix(data = train, label = train_label),
      nrounds = xgb_nrounds,
      params = xgb_params
    )
    new_train_error <- tail(xgb$evaluation_log$train_error, 1)

    if (k >= 2) {
      if (new_train_error > train_error) {
        message("*** train_error increased(", new_train_error, ">", train_error, "). Use the previous model for final classification.")
        break
      }
      if (train_error - new_train_error <= min_improve | new_train_error <= min_error) {
        message("*** train_error is limited(", new_train_error, "). Use the current model for final classification.")
        k <- max_iter
      }
    } else {
      message("*** The first iteration is for rough classification, skipping error checking.")
    }

    train_error <- new_train_error
    model_use <- xgb
    train_use <- train
    train_label_use <- train_label
    XGBoostScore <- predict(model_use, to_predict)
    names(XGBoostScore) <- rownames(to_predict)

    Cell_XGBoostScore <- matrix(XGBoostScore[c(colnames(Cell_counts), colnames(Sim_Cell_counts))], ncol = Cell_downsample_times + 1)
    rownames(Cell_XGBoostScore) <- colnames(Cell_counts)
    Uncertain_XGBoostScore <- matrix(XGBoostScore[colnames(Sim_Uncertain_counts)], ncol = Uncertain_downsample_times)
    rownames(Uncertain_XGBoostScore) <- colnames(Uncertain_counts)
    mu <- abs(cell_score - 0.5)
    stat_list <- lapply(list(Cell_XGBoostScore, Uncertain_XGBoostScore), function(m) {
      pvalue <- apply(m, 1, function(x) {
        if (length(unique(x)) == 1) {
          p <- ifelse(all(abs(x - 0.5) == mu), 1, 0)
        } else {
          p <- t.test(abs(x - 0.5), mu = mu)$p.value
        }
        return(p)
      })
      mean_value <- rowMeans(m)
      return(data.frame(pvalue = pvalue, mean_value = mean_value))
    })
    stat_out <- as.data.frame(Reduce(function(x, y) rbind(x, y), stat_list))
    stat_out[colnames(Uncertain_counts), "FDR"] <- p.adjust(stat_out[colnames(Uncertain_counts), "pvalue"], "BH")

    XGBoostScore <- c(
      stat_out[colnames(Cell_counts), "mean_value"],
      stat_out[colnames(Uncertain_counts), "mean_value"],
      XGBoostScore[colnames(Empty_counts)]
    )
    names(XGBoostScore) <- c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts))

    if (is.null(preCell_mask)) {
      preCell_mask <- train_error > 0.1
    }
    if (isTRUE(preCell_mask)) {
      XGBoostScore[colnames(Cell_counts)] <- 1
    }

    multiscore <- XGBoostScore
    if (isTRUE(Gini_control)) {
      raw_score <- multiscore
      multiscore <- (multiscore + ifelse(multiscore > 0.5, 1, -1) * meta_info[names(multiscore), "CellGiniScore"] + ifelse(multiscore > 0.5, 0, 1)) / 2
      multiscore <- ifelse(raw_score > 0.5, pmin(raw_score, multiscore), pmax(raw_score, multiscore))
    }

    meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "XGBoostScore"] <- XGBoostScore
    meta_info[, "pvalue"] <- 1
    meta_info[, "FDR"] <- 1
    meta_info[rownames(stat_out), "pvalue"] <- stat_out[rownames(stat_out), "pvalue"]
    meta_info[colnames(Uncertain_counts), "FDR"] <- stat_out[colnames(Uncertain_counts), "FDR"]
    meta_info[, paste0("dropSplitScore_iter", j)] <- 0.5
    meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), paste0("dropSplitScore_iter", j)] <- multiscore

    if (j == 1) {
      meta_info[, "dropSplitClass_pre"] <- meta_info[, "dropSplitClass"] <- meta_info[, "preDefinedClass"]
      meta_info[, "dropSplitScore"] <- meta_info[, "dropSplitScore_iter1"]
    } else {
      ## consistent and improved
      meta_info[, "dropSplitClass_pre"] <- meta_info[, "dropSplitClass"]
      cons <- ifelse(meta_info[, "dropSplitScore_iter1"] <= 0.5, 1, -1) * ifelse(meta_info[, paste0("dropSplitScore_iter", j)] <= 0.5, 1, -1)
      improved <- ifelse(abs(meta_info[, paste0("dropSplitScore_iter", j)] - 0.5) > abs(meta_info[, paste0("dropSplitScore_iter", j - 1)] - 0.5), 1, 0)
      meta_info[, "dropSplitScore"] <- ifelse(cons * improved == 1, meta_info[, paste0("dropSplitScore_iter", j)], meta_info[, "dropSplitScore"])
    }

    ## make classification
    meta_info[, "dropSplitClass"] <- ifelse(
      meta_info[, "dropSplitScore"] > cell_score, "Cell", ifelse(meta_info[, "dropSplitScore"] < empty_score, "Empty", "Uncertain")
    )

    ## control FDR for 'Cell' switched from 'Uncertain'
    meta_info[meta_info$preDefinedClass == "Uncertain" & meta_info$dropSplitClass == "Cell" & meta_info$FDR >= cell_score_FDR, "dropSplitClass"] <- "Uncertain"

    ## mask 'Cell' switched from 'Empty'
    if (isTRUE(preEmpty_mask)) {
      meta_info[meta_info$preDefinedClass == "Empty" & meta_info$dropSplitClass == "Cell", "dropSplitClass"] <- "Uncertain"
    }
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))

    message(
      "\n>>> Summary of dropSplit-defined droplet (iteration=", j, ")",
      "\n... #'Cell' Currently defined: ", sum(meta_info$dropSplitClass == "Cell"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Cell"),
      "\n... #'Uncertain' Currently defined: ", sum(meta_info$dropSplitClass == "Uncertain"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Uncertain"),
      "\n... #'Empty' Currently defined: ", sum(meta_info$dropSplitClass == "Empty"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Empty"),
      "\n>>> 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell"),
      "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", "XGBoostScore"]), 3),
      "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", "CellGiniScore"]), 3),
      "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", paste0("dropSplitScore_iter", j)]), 3),
      "\n>>> 'Uncertain' or 'Empty' switch to 'Cell': ", sum(meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell"),
      "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
      "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", "CellGiniScore"]), 3),
      "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", paste0("dropSplitScore_iter", j)]), 3)
    )

    if (do_plot) {
      p <- CellEntropyPlot(meta_info, colorBy = paste0("dropSplitScore_iter", j), splitBy = NULL, cell_stat_by = "dropSplitClass") +
        labs(title = paste0("dropSplitScore(iteration=", j, ")"))
      suppressWarnings(print(p))
    }
    if (j == max_iter) {
      message(">>> Reached maximum number of iterations. Use the last model for final classification.")
    }
  }
  meta_info[, "dropSplitClass_pre"] <- NULL
  if (do_plot) {
    p <- CellEntropyPlot(meta_info, colorBy = "dropSplitScore", splitBy = NULL, cell_stat_by = "dropSplitClass") +
      labs(title = paste0("dropSplitScore(Final)"))
    suppressWarnings(print(p))
  }
  message("\n  ============= Final result =============  ")

  rescure <- which(meta_info[, "preDefinedClass"] %in% c("Uncertain", "Empty") & meta_info[, "dropSplitClass"] == "Cell")
  rescure_score <- meta_info[rescure, "dropSplitScore"]
  er_rate <- train_error / cell_score * 2
  drop <- floor(length(rescure) * er_rate)
  if (drop > 0) {
    drop_index <- rescure[order(rescure_score, decreasing = FALSE)[1:drop]]
    drop_score <- round(mean(meta_info[drop_index, "dropSplitScore"]), digits = 3)
    meta_info[drop_index, "dropSplitScore"] <- 0.5
    meta_info[drop_index, "dropSplitClass"] <- "Uncertain"
  } else {
    drop_score <- NA
  }
  message(
    ">>> Control the false positive rates for 'Cell'",
    "\n... Number of new defined Cell from 'Uncertain' or 'Empty': ", length(rescure),
    "\n... Estimated error rate: ", round(er_rate, digits = 3),
    "\n... Estimated error number: ", drop,
    "\n... Estimated error droplets mean ", "dropSplitScore", ": ", drop_score,
    "\n*** The dropSplitScore and dropSplitClass for these ", drop, " droplets are converted to 0.5 and 'Uncertain'"
  )

  message(
    "\n>>> Summary of final dropSplit-defined droplet",
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
  if (do_plot) {
    p <- CellEntropyPlot(meta_info, colorBy = "dropSplitClass", splitBy = NULL, cell_stat_by = "dropSplitClass") +
      labs(title = paste0("dropSplitClass(Final)"))
    suppressWarnings(print(p))
  }

  importance_matrix <- as.data.frame(xgb.importance(model = model_use))
  result <- list(
    meta_info = DataFrame(meta_info[raw_cell_order, ]),
    train = train_use,
    train_label = train_label_use,
    to_predict = to_predict,
    model = model_use,
    importance_matrix = importance_matrix
  )

  end <- Sys.time()
  message(
    "Elapsed: ", round(difftime(time1 = end, time2 = start, units = "mins"), digits = 3), " mins",
    "\n>>> dropSplit identified ", sum(meta_info$dropSplitClass == "Cell"), " cell-containing droplets."
  )
  return(result)
}

#' Calculate Mean Squared Error for nCount/nFeature Rank.
#' @param meta_info A \code{data.frame} or \code{DataFrame}. Must contain the 'nCount' and 'nFeature' columns.
#' @param smooth_num Number of times to smooth(take a mean value within a window length \code{smooth_window}) the squared error. Default is 5.
#' @param smooth_window Window length used to smooth the squared error. Default is 100.
#' @param find_cell_count Whether to find RankMSE valleys ('Cell' count cutoff). Default is FALSE.
#' @param tolerance A value indicated the tolerance when finding RankMSE valleys. A value greater than 1 indicates relaxed and will find more valleys; lower than 1 indicates strict and will find less valleys. Default is 0.5.
#'
#' @return A list include \code{meta_info} and \code{cell_rank_count}
RankMSE <- function(meta_info, smooth_num = 5, smooth_window = 100, find_cell_count = FALSE, tolerance = 0.5) {
  meta_info <- as.data.frame(meta_info)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  meta_info$nFeature_rank <- rank(-(meta_info$nFeature))

  if (min(meta_info$nCount) <= 0) {
    stop("Find 'nCount'<=0. Stop run.")
  }
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]

  meta_info$RankSE <- (meta_info$nCount_rank - meta_info$nFeature_rank)^2
  x0 <- meta_info$RankSE
  x0_nCount <- meta_info$nCount
  x0_cell <- rownames(meta_info)
  x1 <- rep(x0[-length(x0)], -diff(x0_nCount))
  x1_nCount <- rep(x0_nCount[-length(x0_nCount)], -diff(x0_nCount))
  x1_cell <- rep(x0_cell[-length(x0_cell)], -diff(x0_nCount))
  x1_cell <- paste0("Sim-", seq_len(length(x1_cell)), "@Raw-", x1_cell)
  filled_df <- rbind(
    data.frame(RankMSE = x0, nCount = x0_nCount, cell = x0_cell),
    data.frame(RankMSE = x1, nCount = x1_nCount, cell = x1_cell)
  )
  filled_df <- filled_df[order(filled_df[, "nCount"], decreasing = TRUE), ]

  for (t in seq_len(smooth_num)) {
    filled_df[, "RankMSE"] <- runMean(filled_df[, "RankMSE"], n = smooth_window)
    filled_df[, "RankMSE"][is.na(filled_df[, "RankMSE"])] <- na.omit(filled_df[, "RankMSE"])[1]
    filled_df[, "RankMSE"][filled_df[, "RankMSE"] < 0] <- 0
  }
  filled_df[, "logRankMSE"] <- log(filled_df[, "RankMSE"])

  if (isTRUE(find_cell_count)) {
    nCount_inflection <- tail(find_inflection(meta_info$nCount)$index, 1)
    nFeature_inflection <- tail(find_inflection(meta_info$nFeature)$index, 1)
    inflection <- max(nCount_inflection, nFeature_inflection)
    inflection_left <- round(inflection - inflection * 0.3)
    inflection_right <- round(inflection + inflection * 0.3)
    inflection_left <- head(which(filled_df[, "nCount"] <= meta_info$nCount[inflection_left]), 1)
    inflection_right <- tail(which(filled_df[, "nCount"] >= meta_info$nCount[inflection_right]), 1)
    pks <- inflection_left + find_peaks(-filled_df$RankMSE[inflection_left:inflection_right],
      left_shoulder = inflection_left * 0.05,
      right_shoulder = inflection_right
    ) - 1
    # qplot((1:length(filled_df$RankMSE)), log(filled_df$RankMSE)) + xlim(1, inflection_right*1.2) + geom_vline(xintercept = c(pks))

    iter <- 1
    cell_max <- quantile(filled_df[1:(pks[1] - 1), "logRankMSE"], 0.99)
    while (length(pks) > 1) {
      message("... ", length(pks), " peaks found. Perform automatic selection(iter=", iter, ")...")
      if (all(filled_df[pks, "logRankSE"] < 0)) {
        pks <- pks[length(pks)]
      }
      if (length(pks) > 1) {
        pks_diff_MSE <- diff(filled_df[pks, "logRankMSE"])
        pks_max_MSE <- sapply(1:(length(pks) - 1), function(i) {
          quantile(filled_df[pks[i]:pks[i + 1], "logRankMSE"], 0.99) - filled_df[pks[i + 1], "logRankMSE"]
        })
        j <- which(pks_diff_MSE <= pks_max_MSE * tolerance | pks_max_MSE >= 0.5 * cell_max) + 1
        if (length(j) == 0) {
          j <- 1
        }
        pks <- pks[j]
        iter <- iter + 1
      }
      if (length(pks) == 1) {
        message("... Automatic selection finished.")
      } else {
        message("... ", length(pks), " apples left. Go to the next loop.")
      }
    }
    cell_count <- filled_df[pks, "nCount"]
  } else {
    cell_count <- NULL
  }

  raw_df <- unique(filled_df[, c("cell", "RankMSE")])
  rownames(raw_df) <- raw_df$cell
  meta_info$RankMSE <- raw_df[rownames(meta_info), "RankMSE"]

  return(list(meta_info = meta_info, cell_count = cell_count))
}

CalldropSplit <- function(counts, ...) {
  result <- dropSplit(counts, ...)
  return(result)
}
