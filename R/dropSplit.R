#' Automatically identify cell-containing and empty droplets for droplet-based scRNAseq data using dropSplit.
#'
#' @description dropSplit is designed to identify true cells from droplet-based scRNAseq data.
#' It consists of four main steps:
#' \enumerate{
#' \item Pre-define droplets as 'Cell', 'Uncertain', 'Empty' and 'Discarded' droplets according to the RankMSE curve.
#' \item Simulate 'Cell' and 'Uncertain' droplets under a depth of 'Empty' used for model construction and prediction.
#' \item Iteratively buid the model, classify 'Uncertain' droplets, and update the training set use newly predicted 'Empty'.
#' \item Make classification with the optimal model.
#' }
#' dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini which can help identification.
#' In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification.
#' It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent droplets and rows represent features.
#' @param do_plot Whether to plot during the cellcalling. Default is \code{TRUE}.
#' @param Cell_score A cutoff value of \code{dropSplitScore} to determine if a droplet is cell-containing. Range between 0.5 and 1. Default is 0.9.
#' @param Empty_score A cutoff value of \code{dropSplitScore} to determine if a droplet is empty. Range between 0 and 0.5. Default is 0.2.
#' @param downsample_times Number of repetitions of downsampling for 'Cell' and 'Uncertain' droplets. If \code{NULL}, will be determined by \code{CE_ratio}. Default is \code{NULL}.
#' @param CE_ratio Ratio value between down-sampled 'Cells' and 'Empty' droplets. The actual value will be slightly higher than the set. Default is 1.
#' @inheritParams RankMSE
#' @param Cell_min_nCount Minimum nCount for 'Cell' droplets. Default is 500.
#' @param Empty_min_nCount Minimum nCount for 'Empty' droplets. Default is 20.
#' @param Empty_max_num Number of pre-defined 'Empty' droplets. Default is 50000.
#' @param max_iter An integer specifying the number of iterations to use to rebuild the model with new defined droplets. Default is 5.
#' @param Gini_control Whether to control cell quality by CellGini. Default is \code{TRUE}.
#' @param Gini_threshold A cutoff of the CellGini metric. The higher, the more conservative and will get a lower number of cells. Default is automatic.
#' @param Cell_rank,Uncertain_rank,Empty_rank Custom Rank value to mark the droplets as Cell, Uncertain and Empty labels for the data to be trained. Default is automatic. But useful when the default value is considered to be wrong from the RankMSE plot.
#' @param preCell_mask logical; Whether to mask pre-defined 'Cell' droplets when prediction. If \code{TRUE}, XGBoostScore for all droplets pre-defined as 'Cell' will be set to 1; Default is \code{FALSE}.
#' @param preEmpty_mask logical; Whether to mask pre-defined 'Empty' droplets when prediction. There is a little different with parameter \code{preCell_mask}. If \code{TRUE}, XGBoostScore will not change, but the final classification will not be 'Cell' in any case. Default is \code{TRUE}.
#' @param FDR FDR cutoff for droplets that predicted as 'Cell' or 'Empty' from pre-defined 'Uncertain'. Note, statistic tests and the FDR control only performed on the difference between averaged \code{XGBoostScore} and 0.5. Default is 0.05.
#' @param xgb_params The \code{list} of XGBoost parameters.
#' @param modelOpt Whether to optimize the model using \code{\link{xgbOptimization}}. Will take long time for large datasets. If \code{TRUE}, will overwrite the parameters list in \code{xgb_params}. The following parameters are only used in \code{\link{xgbOptimization}}.
#' @inheritParams xgbOptimization
#' @param verbose If 0, xgboost will stay silent. If 1, it will print information about performance. If 2, some additional information will be printed out. Note that setting verbose > 0 automatically engages the cb.print.evaluation(period=1) callback function.
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
#' @importFrom inflection uik
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
dropSplit <- function(counts, do_plot = TRUE, Cell_score = 0.9, Empty_score = 0.2,
                      downsample_times = NULL, CE_ratio = 1,
                      fill_RankMSE = FALSE, smooth_num = 2, smooth_window = 100, tolerance = 0.2,
                      Cell_rank = NULL, Uncertain_rank = NULL, Empty_rank = NULL,
                      Cell_min_nCount = 500, Empty_min_nCount = 20, Empty_max_num = 50000,
                      Gini_control = TRUE, Gini_threshold = NULL, max_iter = 5,
                      preCell_mask = FALSE, preEmpty_mask = TRUE, FDR = 0.05,
                      xgb_params = NULL, xgb_nrounds = 20, xgb_thread = 8, xgb_early_stopping_rounds = 3,
                      modelOpt = FALSE, verbose = 1, seed = 0, ...) {
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
  if (!is.null(Cell_score)) {
    if (Cell_score < 0.5 | Cell_score > 1) {
      stop("'Cell_score' must be between 0.5 and 1.")
    }
  }
  if (!is.null(Gini_threshold)) {
    if (Gini_threshold < 0 | Gini_threshold > 1) {
      stop("'Gini_threshold' must be between 0 and 1.")
    }
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

  out <- RankMSE(
    meta_info = meta_info, fill_RankMSE = fill_RankMSE, smooth_num = smooth_num, smooth_window = smooth_window,
    find_rank = TRUE, tolerance = tolerance, Empty_min_nCount = Empty_min_nCount
  )
  meta_info <- out$meta_info
  inflection <- out$inflection
  Cell_rank <- out$Cell_rank
  Uncertain_rank <- out$Uncertain_rank
  Empty_rank <- out$Empty_rank

  message(
    "... The inflection point is detected at ", inflection,
    "\n... Cell_rank=", Cell_rank,
    "\n... Uncertain_rank=", Uncertain_rank,
    "\n... Empty_rank=", Empty_rank
  )
  Cell_count <- meta_info$nCount[Cell_rank]
  Empty_count <- meta_info$nCount[Empty_rank]
  Uncertain_count <- meta_info$nCount[Uncertain_rank]

  Cell_counts <- counts[, meta_info$nCount >= Cell_count]
  Empty_counts <- counts[, meta_info$nCount < Uncertain_count & meta_info$nCount >= Empty_count]
  Uncertain_counts <- counts[, meta_info$nCount < Cell_count & meta_info$nCount >= Uncertain_count]

  if (do_plot) {
    p <- RankMSEPlot(meta_info, colorBy = "nFeature", splitBy = NULL, cell_stat_by = NULL) +
      geom_vline(
        xintercept = c(inflection),
        color = c("black")
      ) +
      geom_vline(
        xintercept = c(Cell_rank, Uncertain_rank, Empty_rank),
        color = c("red3", "forestgreen", "steelblue"), linetype = 2
      )
    suppressWarnings(print(p))
  }

  if (ncol(Uncertain_counts) == 0) {
    stop("No 'Uncertain' droplets detected. Please check the RankMSE curve with the pre-defined droplet cutoff. You may set custom cutoff values in the parameters manually.")
  }
  if (ncol(Empty_counts) == 0) {
    stop("No 'Empty' droplets detected. Please check the RankMSE curve with the pre-defined droplet cutoff. You may set custom cutoff values in the parameters manually.")
  }
  if (ncol(Empty_counts) < ncol(Cell_counts)) {
    stop("Pre-defined 'Empty' droplets is fewer than 'Cell'. You may set custom rank values in the parameters manually or set a lower value for 'Empty_min_nCount'.")
  }
  if (ncol(Empty_counts) > Empty_max_num) {
    warning("The number of 'Empty' droplets exceeds the specified. Converting the excess 'Empty' droplets to 'Uncertain' droplets.", immediate. = TRUE)
    Uncertain_counts <- cbind(Uncertain_counts, Empty_counts[, 1:(ncol(Empty_counts) - Empty_max_num)])
    Empty_counts <- Empty_counts[, (ncol(Empty_counts) - Empty_max_num + 1):ncol(Empty_counts)]
  }

  message(">>> Calculate CellGiniScore for the droplets...")
  primary_counts <- cbind(Cell_counts, Uncertain_counts, Empty_counts)
  final_Gini <- CellGini(primary_counts, normalize = TRUE)
  if (is.null(Gini_threshold)) {
    Gini_threshold <- max(
      quantile(final_Gini[colnames(Cell_counts)], 0.25) - diff(quantile(final_Gini[colnames(Cell_counts)], c(0.25, 0.75))) * 3,
      min(final_Gini[colnames(Cell_counts)])
    )
  }
  message("*** Gini_threshold was set to: ", round(Gini_threshold, 3))
  minGini <- max(Gini_threshold - 0.05, min(final_Gini))
  final_CellGiniScore <- (final_Gini - minGini) / (Gini_threshold - minGini)
  final_CellGiniScore[final_CellGiniScore < 0] <- 0
  final_CellGiniScore[final_CellGiniScore > 1] <- 1
  names(final_CellGiniScore) <- names(final_Gini)
  dat_Features <- cbind(
    CellGini = final_Gini,
    CellGiniScore = final_CellGiniScore
  )
  rownames(dat_Features) <- colnames(primary_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_Features)
  ] <- dat_Features[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]

  cell_Gini <- final_Gini[colnames(Cell_counts)]
  cell_drop <- names(cell_Gini)[cell_Gini < Gini_threshold]
  cell_use <- names(cell_Gini)[cell_Gini >= Gini_threshold]
  if (length(cell_drop) > 0) {
    message(
      "*** There are ", length(cell_drop), " high RNA content droplets with low 'CellGini' values than the Gini_threshold:", round(Gini_threshold, 3), "\n",
      "*** These ", length(cell_drop), " droplets will pre-defined as 'Uncertain'"
    )
    Drop_counts <- Cell_counts[, cell_drop]
    if (is.numeric(Drop_counts)) {
      Drop_counts <- as.matrix(Drop_counts)
      colnames(Drop_counts) <- cell_drop
    }
    Cell_counts <- Cell_counts[, cell_use]
    Uncertain_counts <- cbind(as(Drop_counts, "dgCMatrix"), Uncertain_counts)
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
    p <- RankMSEPlot(meta_info, colorBy = "preDefinedClass", splitBy = NULL, cell_stat_by = "preDefinedClass") + geom_vline(xintercept = inflection)
    suppressWarnings(print(p))
  }

  message(
    ">>> Downsample pre-defined 'Cell' droplets to a depth similar to the 'Empty' for training",
    "\n... nCount in 'Cell' droplets: Min=", min(Cell_nCount), " Median=", median(Cell_nCount), " Max=", max(Cell_nCount),
    "\n... nCount in 'Empty' droplets: Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount)
  )
  if (!is.null(CE_ratio)) {
    downsample_times <- ceiling(ncol(Empty_counts) / ncol(Cell_counts) * CE_ratio) - 1
  }
  Sim_Cell_counts <- Cell_counts[, rep(1:ncol(Cell_counts), downsample_times)]
  colnames(Sim_Cell_counts) <- paste0(rep(paste0("Sim", 1:downsample_times), each = ncol(Cell_counts)), "-", colnames(Sim_Cell_counts))
  Sim_Cell_nCount_assign <- sample(Empty_nCount, ncol(Sim_Cell_counts), replace = TRUE)
  Sim_Cell_counts <- downsampleMatrix(x = Sim_Cell_counts, prop = Sim_Cell_nCount_assign / Matrix::colSums(Sim_Cell_counts), bycol = TRUE)
  Sim_Cell_nCount <- Matrix::colSums(Sim_Cell_counts)
  CE_ratio <- (ncol(Cell_counts) + ncol(Sim_Cell_counts)) / ncol(Empty_counts)
  message(">>> 'Cell/Empty' ratio is adjusted to ", round(CE_ratio, 3), "(downsample_times=", downsample_times, ")")

  message(
    ">>> Downsample pre-defined 'Uncertain' droplets to a depth similar to the 'Empty' for prediction",
    "\n... nCount in 'Uncertain' droplets: Min=", min(Uncertain_nCount), " Median=", median(Uncertain_nCount), " Max=", max(Uncertain_nCount),
    "\n... nCount in 'Empty' droplets: Min=", min(Empty_nCount), " Median=", median(Empty_nCount), " Max=", max(Empty_nCount)
  )
  Sim_Uncertain_counts <- Uncertain_counts[, rep(1:ncol(Uncertain_counts), downsample_times)]
  colnames(Sim_Uncertain_counts) <- paste0(rep(paste0("Sim", 1:downsample_times), each = ncol(Uncertain_counts)), "-", colnames(Sim_Uncertain_counts))
  Sim_Uncertain_nCount_assign <- sample(Empty_nCount, ncol(Sim_Uncertain_counts), replace = TRUE)
  Sim_Uncertain_counts <- downsampleMatrix(x = Sim_Uncertain_counts, prop = Sim_Uncertain_nCount_assign / Matrix::colSums(Sim_Uncertain_counts), bycol = TRUE)
  Sim_Uncertain_nCount <- Matrix::colSums(Sim_Uncertain_counts)

  message(">>> Calculate various metrics for the droplets to be trained...")
  comb_counts <- cbind(Cell_counts, Sim_Cell_counts, Uncertain_counts, Sim_Uncertain_counts, Empty_counts)
  MTgene <- grep(x = rownames(comb_counts), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T, value = TRUE)
  RPgene <- grep(x = rownames(comb_counts), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T, value = TRUE)
  comb_MTprop <- Matrix::colSums(comb_counts[MTgene, ]) / Matrix::colSums(comb_counts)
  comb_RPprop <- Matrix::colSums(comb_counts[RPgene, ]) / Matrix::colSums(comb_counts)
  comb_CellEntropy <- CellEntropy(comb_counts)
  comb_maxCellEntropy <- maxCellEntropy(comb_counts)
  comb_CellEfficiency <- comb_CellEntropy / comb_maxCellEntropy
  comb_CellEfficiency[is.na(comb_CellEfficiency)] <- 1
  comb_CellRedundancy <- 1 - comb_CellEfficiency
  comb_CellRedundancy[comb_CellRedundancy < 0] <- 0
  nCount_per_Feature <- Matrix::colSums(comb_counts) / Matrix::colSums(comb_counts > 0)
  dat_Features <- cbind(
    MTprop = comb_MTprop,
    RPprop = comb_RPprop,
    CellEntropy = comb_CellEntropy,
    CellRedundancy = comb_CellRedundancy,
    nCount_per_Feature = nCount_per_Feature
  )
  rownames(dat_Features) <- colnames(comb_counts)
  meta_info[
    c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)),
    colnames(dat_Features)
  ] <- dat_Features[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), ]

  ### create bins for 'Uncertain'
  uncertain_maxCount <- min(Cell_nCount) - 1
  uncertain_minCount <- min(Uncertain_nCount)
  breaks <- c(0, 10^seq(log10(uncertain_minCount), log10(uncertain_maxCount), length.out = max_iter), max(meta_info$nCount))
  meta_info[, "nCount_Bin"] <- as.numeric(cut(meta_info$nCount, breaks = breaks))
  meta_info[colnames(Empty_counts), "nCount_Bin"] <- 0
  meta_info[colnames(Cell_counts), "nCount_Bin"] <- 1

  message(">>> Merge new features into train data...")
  norm_counts <- Matrix::t(Matrix::t(comb_counts) / Matrix::colSums(comb_counts))
  dat <- cbind(Matrix::t(norm_counts), dat_Features[, colnames(dat_Features)])
  ini_train <- dat[c(colnames(Cell_counts), colnames(Sim_Cell_counts), colnames(Empty_counts)), ]
  ini_train_label <- c(rep(1, ncol(Cell_counts) + ncol(Sim_Cell_counts)), rep(0, ncol(Empty_counts)))
  ini_to_predict <- dat[c(
    colnames(Cell_counts), colnames(Sim_Cell_counts),
    colnames(Uncertain_counts), colnames(Sim_Uncertain_counts),
    colnames(Empty_counts)
  ), ]
  Sim_Cell_counts_rawname <- gsub(x = colnames(Sim_Cell_counts), pattern = "Sim\\d+-", replacement = "")
  Sim_Uncertain_counts_rawname <- gsub(x = colnames(Sim_Uncertain_counts), pattern = "Sim\\d+-", replacement = "")

  xgb_base_params <- list(
    eta = 0.3,
    gamma = 1,
    max_depth = 20,
    min_child_weight = 7,
    subsample = 0.7,
    max_delta_step = 1,
    alpha = 5,
    lambda = 10,
    eval_metric = "error",
    eval_metric = "aucpr",
    objective = "binary:logistic",
    nthread = xgb_thread
  )
  if (!is.null(xgb_params)) {
    for (nm in names(xgb_params)) {
      xgb_base_params[[nm]] <- xgb_params[[nm]]
    }
  }
  xgb_params <- xgb_base_params

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
  message(">>> Construct the XGBoost model with pre-defined classification...")

  train_error <- 1
  k <- 0
  j <- 0
  while (k < max_iter) {
    k <- k + 1
    j <- j + 1
    message("\n  ============= Iteration: ", k, " =============  ")
    if (k == 1) {
      train <- ini_train
      train_label <- ini_train_label
      to_predict <- ini_to_predict
      empty_update <- colnames(Empty_counts)
      cell_update <- c(colnames(Cell_counts), colnames(Sim_Cell_counts))
      highlight <- colnames(Cell_counts)
    } else {
      ### empty
      if (k == 2) {
        empty_current <- rownames(meta_info)[meta_info$dropSplitClass == "Empty" & meta_info$nCount_Bin %in% c(0, k)]
      } else {
        empty_current <- c(empty_update, rownames(meta_info)[meta_info$dropSplitClass == "Empty" & meta_info$nCount_Bin == k])
      }
      new_empty <- colnames(Sim_Uncertain_counts)[Sim_Uncertain_counts_rawname %in% empty_current]
      if (length(new_empty) > ceiling(Empty_max_num / max_iter)) {
        new_empty <- sample(new_empty, ceiling(Empty_max_num / max_iter))
      }
      empty_remove_num <- length(empty_update[empty_update %in% rownames(meta_info)]) - length(empty_update[empty_update %in% empty_current & empty_update %in% rownames(meta_info)])
      raw_empty <- empty_update[empty_update %in% empty_current]
      empty_update <- unique(c(new_empty, raw_empty))

      ### cell
      if (k == 2) {
        cell_current <- rownames(meta_info)[meta_info$dropSplitClass == "Cell" & meta_info$nCount_Bin %in% c(1, max_iter)]
      } else {
        cell_current <- c(cell_update, rownames(meta_info)[meta_info$dropSplitClass == "Cell" & meta_info$nCount_Bin %in% c(1, max_iter - k + 2)])
      }
      new_cell <- colnames(Sim_Uncertain_counts)[Sim_Uncertain_counts_rawname %in% cell_current]
      if (length(new_cell) > ceiling(Empty_max_num / max_iter)) {
        new_cell <- sample(new_cell, ceiling(Empty_max_num / max_iter))
      }
      cell_remove_num <- length(cell_update[cell_update %in% rownames(meta_info)]) - length(cell_update[cell_update %in% cell_current & cell_update %in% rownames(meta_info)])
      raw_cell <- cell_update[cell_update %in% cell_current]
      cell_update <- unique(c(new_cell, raw_cell))
      if (length(cell_update) / length(empty_update) < CE_ratio) {
        nsample <- ceiling(length(empty_update) * CE_ratio - length(cell_update))
        new_cell <- sample(cell_update, nsample, replace = TRUE)
        cell_update <- c(cell_update, new_cell)
      }

      if (length(empty_update) < ncol(Cell_counts)) {
        message("*** Number of 'Empty' droplets is too small(", length(empty_update), "). Use the previous model for final classification.")
        break
      }
      if (length(new_empty) == 0 & length(new_cell) == 0 & empty_remove_num == 0 & cell_remove_num == 0) {
        message("*** No change in the training data. Use the previous model for final classification.")
        break
      }

      train <- dat[c(cell_update, empty_update), ]
      train_label <- c(rep(1, length(cell_update)), rep(0, length(empty_update)))

      highlight <- intersect(
        cell_current,
        rownames(meta_info)[meta_info[, "preDefinedClass"] == "Uncertain"]
      )
      message(
        "... Number of filtered 'Cell' droplets defined in the previous round: ", cell_remove_num,
        "\n... Number of filtered 'Empty' droplets defined in the previous round: ", empty_remove_num,
        "\n... Number of the current training 'Cell' droplets from 'Uncertain': ", length(intersect(
          cell_current,
          rownames(meta_info)[meta_info[, "preDefinedClass"] == "Uncertain"]
        )),
        "\n... Number of the current training 'Empty' droplets from 'Uncertain': ", length(intersect(
          empty_current,
          rownames(meta_info)[meta_info[, "preDefinedClass"] == "Uncertain"]
        ))
      )
    }

    message(">>> Train data: 'Cell'=", sum(train_label == 1), "; 'Empty'=", sum(train_label == 0), "; 'C/E' Ratio=", round(sum(train_label == 1) / sum(train_label == 0), 3))
    xgb <- xgboost(
      data = xgb.DMatrix(data = train, label = train_label),
      nrounds = xgb_nrounds,
      params = xgb_params,
      early_stopping_rounds = xgb_early_stopping_rounds,
      verbose = verbose
    )
    new_train_error <- tail(xgb$evaluation_log$train_error, 1)
    new_train_aucpr <- tail(xgb$evaluation_log$train_aucpr, 1)

    train_error <- new_train_error
    train_aucpr <- new_train_aucpr
    model_use <- xgb
    train_use <- train
    train_label_use <- train_label
    message(">>> Predict the type of droplets...")
    XGBoostScore <- predict(model_use, to_predict)
    names(XGBoostScore) <- rownames(to_predict)

    Cell_XGBoostScore <- matrix(XGBoostScore[c(colnames(Sim_Cell_counts))], ncol = downsample_times)
    rownames(Cell_XGBoostScore) <- colnames(Cell_counts)
    Uncertain_XGBoostScore <- matrix(XGBoostScore[c(colnames(Sim_Uncertain_counts))], ncol = downsample_times)
    rownames(Uncertain_XGBoostScore) <- colnames(Uncertain_counts)
    # mu <- abs(Cell_score - 0.5)
    mu <- 0
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

    xgCellScore <- ifelse(XGBoostScore[colnames(Cell_counts)] > 0.5,
      pmax(stat_out[colnames(Cell_counts), "mean_value"], XGBoostScore[colnames(Cell_counts)]),
      pmin(stat_out[colnames(Cell_counts), "mean_value"], XGBoostScore[colnames(Cell_counts)])
    )
    xgUncertainScore <- stat_out[colnames(Uncertain_counts), "mean_value"]
    XGBoostScore <- c(
      xgCellScore,
      xgUncertainScore,
      XGBoostScore[colnames(Empty_counts)]
    )
    names(XGBoostScore) <- c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts))

    if (isTRUE(preCell_mask)) {
      XGBoostScore[colnames(Cell_counts)] <- 1
    }

    multiscore <- XGBoostScore
    if (isTRUE(Gini_control)) {
      raw_score <- multiscore
      multiscore <- (multiscore + ifelse(multiscore > 0.5, 1, -1) * meta_info[names(multiscore), "CellGiniScore"] + ifelse(multiscore > 0.5, 0, 1)) / 2
      multiscore <- ifelse(raw_score > 0.5, pmin(raw_score, multiscore), pmax(raw_score, multiscore))
      multiscore <- ifelse(raw_score > 0.5, pmax(0.5, multiscore), pmin(0.5, multiscore))
    }

    meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), "XGBoostScore"] <- XGBoostScore
    meta_info[, "pvalue"] <- 1
    meta_info[, "FDR"] <- 1
    meta_info[rownames(stat_out), "pvalue"] <- stat_out[rownames(stat_out), "pvalue"]
    meta_info[colnames(Uncertain_counts), "FDR"] <- stat_out[colnames(Uncertain_counts), "FDR"]
    meta_info[, paste0("dropSplitScore_iter", j)] <- 0.5
    meta_info[c(colnames(Cell_counts), colnames(Uncertain_counts), colnames(Empty_counts)), paste0("dropSplitScore_iter", j)] <- multiscore

    ## make classification
    meta_info[, paste0("dropSplitClass_iter", j)] <- ifelse(
      meta_info[, paste0("dropSplitScore_iter", j)] >= Cell_score, "Cell", ifelse(meta_info[, paste0("dropSplitScore_iter", j)] <= Empty_score, "Empty", "Uncertain")
    )
    meta_info[meta_info[, "preDefinedClass"] == "Discarded", paste0("dropSplitClass_iter", j)] <- "Discarded"

    ## control FDR for 'Cell' switched from 'Uncertain'
    FDR_filter1 <- rownames(meta_info)[meta_info[, "preDefinedClass"] == "Uncertain" & meta_info[, paste0("dropSplitClass_iter", j)] == "Cell" & meta_info$FDR >= FDR]
    message(">>> Filter out ", length(FDR_filter1), " potential 'Cell' droplets by FDR control", "\n... mean FDR: ", round(mean(meta_info[FDR_filter1, "FDR"]), 3))
    meta_info[FDR_filter1, paste0("dropSplitClass_iter", j)] <- "Uncertain"
    ## control FDR for 'Empty' switched from 'Uncertain'
    FDR_filter2 <- rownames(meta_info)[meta_info[, "preDefinedClass"] == "Uncertain" & meta_info[, paste0("dropSplitClass_iter", j)] == "Empty" & meta_info$FDR >= FDR]
    message(">>> Filter out ", length(FDR_filter2), " potential 'Empty' droplets by FDR control", "\n... mean FDR: ", round(mean(meta_info[FDR_filter2, "FDR"]), 3))
    meta_info[FDR_filter2, paste0("dropSplitClass_iter", j)] <- "Uncertain"

    ## mask 'Cell' switched from 'Empty'
    if (isTRUE(preEmpty_mask)) {
      meta_info[meta_info[, "preDefinedClass"] == "Empty" & meta_info[, paste0("dropSplitClass_iter", j)] == "Cell", paste0("dropSplitClass_iter", j)] <- "Uncertain"
    }
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"], levels = c("Cell", "Uncertain", "Empty", "Discarded"))
    meta_info[, paste0("dropSplitClass_iter", j)] <- factor(meta_info[, paste0("dropSplitClass_iter", j)], levels = c("Cell", "Uncertain", "Empty", "Discarded"))

    ## mask 'Cell' with nCount<Cell_min_nCount
    meta_info[meta_info[, paste0("dropSplitClass_iter", j)] == "Cell" & meta_info[, "nCount"] < Cell_min_nCount, paste0("dropSplitClass_iter", j)] <- "Uncertain"

    meta_info[, "dropSplitClass"] <- meta_info[, paste0("dropSplitClass_iter", j)]
    meta_info[, "dropSplitScore"] <- meta_info[, paste0("dropSplitScore_iter", j)]
    if (j == 1) {
      meta_info[, "dropSplitClass_pre"] <- meta_info[, "preDefinedClass"]
    } else {
      meta_info[, "dropSplitClass_pre"] <- meta_info[, paste0("dropSplitClass_iter", j - 1)]
    }

    message(
      "\n>>> Summary of dropSplit-defined droplet (iteration=", j, ")",
      "\n... #'Cell' Currently defined: ", sum(meta_info$dropSplitClass == "Cell"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Cell"),
      "\n... ... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
      "\n... ... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Cell", "CellGiniScore"]), 3),
      "\n... ... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3),
      "\n... #'Uncertain' Currently defined: ", sum(meta_info$dropSplitClass == "Uncertain"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Uncertain"),
      "\n... ... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Uncertain", "XGBoostScore"]), 3),
      "\n... ... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Uncertain", "CellGiniScore"]), 3),
      "\n... ... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Uncertain", "dropSplitScore"]), 3),
      "\n... #'Empty' Currently defined: ", sum(meta_info$dropSplitClass == "Empty"), "  Previously defined: ", sum(meta_info$dropSplitClass_pre == "Empty"),
      "\n... ... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Empty", "XGBoostScore"]), 3),
      "\n... ... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Empty", "CellGiniScore"]), 3),
      "\n... ... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass == "Empty", "dropSplitScore"]), 3),
      "\n>>> 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell"),
      "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", "XGBoostScore"]), 3),
      "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", "CellGiniScore"]), 3),
      "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre == "Cell" & meta_info$dropSplitClass != "Cell", "dropSplitScore"]), 3),
      "\n>>> 'Uncertain' or 'Empty' switch to 'Cell': ", sum(meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell"),
      "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
      "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", "CellGiniScore"]), 3),
      "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info$dropSplitClass_pre != "Cell" & meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3)
    )

    if (do_plot) {
      p <- CellEntropyPlot(meta_info, colorBy = paste0("dropSplitScore_iter", j), splitBy = NULL, cell_stat_by = "dropSplitClass", highlight = highlight) +
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
  er_rate <- train_error * (1 - Cell_score) * 2
  drop <- round(length(rescure) * er_rate)
  if (drop > 0) {
    drop_index <- rescure[order(rescure_score, decreasing = FALSE)[1:drop]]
    drop_score <- round(mean(meta_info[drop_index, "dropSplitScore"]), digits = 3)
    meta_info[drop_index, "dropSplitClass"] <- "Uncertain"
  } else {
    drop_score <- NA
  }
  message(
    ">>> Estimated the classification error for 'Cell'",
    "\n... Number of new defined Cell from 'Uncertain' or 'Empty': ", length(rescure),
    "\n... Estimated error rate: ", round(er_rate, digits = 3),
    "\n... Estimated error number: ", drop,
    "\n... Estimated error droplets mean ", "dropSplitScore", ": ", drop_score,
    "\n*** The dropSplitClass for these ", drop, " droplets are set to 'Uncertain'"
  )

  message(
    "\n>>> Summary of final dropSplit-defined droplet",
    "\n... Number of 'Cell': ", sum(meta_info$dropSplitClass == "Cell"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Cell") > 1, min(meta_info[meta_info$dropSplitClass == "Cell", "nCount"]), 1),
    "\n... Number of 'Uncertain': ", sum(meta_info$dropSplitClass == "Uncertain"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Uncertain") > 1, min(meta_info[meta_info$dropSplitClass == "Uncertain", "nCount"]), 1),
    "\n... Number of 'Empty': ", sum(meta_info$dropSplitClass == "Empty"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Empty") > 1, min(meta_info[meta_info$dropSplitClass == "Empty", "nCount"]), 1),
    "\n... Number of 'Discarded': ", sum(meta_info$dropSplitClass == "Discarded"), "  Minimum nCounts: ", ifelse(sum(meta_info$dropSplitClass == "Discarded") > 1, min(meta_info[meta_info$dropSplitClass == "Discarded", "nCount"]), 1),
    "\n>>> Pre-defined as 'Cell' switch to 'Empty' or 'Uncertain': ", sum(meta_info[, "preDefinedClass"] == "Cell" & meta_info$dropSplitClass != "Cell"),
    "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] == "Cell" & meta_info$dropSplitClass != "Cell", "XGBoostScore"]), 3),
    "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] == "Cell" & meta_info$dropSplitClass != "Cell", "CellGiniScore"]), 3),
    "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] == "Cell" & meta_info$dropSplitClass != "Cell", "dropSplitScore"]), 3),
    "\n>>> Pre-defined as 'Uncertain' or 'Empty' switch to 'Cell': ", sum(meta_info[, "preDefinedClass"] != "Cell" & meta_info$dropSplitClass == "Cell"),
    "\n... Mean XGBoostScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] != "Cell" & meta_info$dropSplitClass == "Cell", "XGBoostScore"]), 3),
    "\n... Mean CellGiniScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] != "Cell" & meta_info$dropSplitClass == "Cell", "CellGiniScore"]), 3),
    "\n... Mean dropSplitScore: ", round(mean(meta_info[meta_info[, "preDefinedClass"] != "Cell" & meta_info$dropSplitClass == "Cell", "dropSplitScore"]), 3)
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
#' @param fill_RankMSE Whether to fill the RankMSE by nCount. Default is \code{TRUE}.
#' @param smooth_num Number of times to smooth(take a mean value within a window length \code{smooth_window}) the squared error. Default is 2.
#' @param smooth_window Window length used to smooth the squared error. Default is 100.
#' @param find_rank Whether to find the 'Cell' RankMSE valley, the 'Uncertain' RankMSE peak and the 'Empty' RankMSE valley. Default is FALSE.
#' @param tolerance A value indicated the tolerance when finding RankMSE valleys. A value greater than 1 indicates relaxed and will find more valleys; lower than 1 indicates strict and will find less valleys. Default is 0.2.
#' @param Empty_min_nCount Minimum nCount for 'Empty' droplets. Default is 20.
#'
#' @return A list include \code{meta_info} and \code{cell_rank_count}
#' @importFrom inflection uik
#' @importFrom TTR runMean
RankMSE <- function(meta_info, fill_RankMSE = FALSE, smooth_num = 2, smooth_window = 100,
                    find_rank = FALSE, tolerance = 0.2, Empty_min_nCount = 20) {
  meta_info <- as.data.frame(meta_info)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  meta_info$nFeature_rank <- rank(-(meta_info$nFeature))

  if (min(meta_info$nCount) <= 0) {
    stop("Find 'nCount'<=0. Stop run.")
  }
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  nCount_inflection <- find_inflection(meta_info$nCount)$index[1]
  nFeature_inflection <- find_inflection(meta_info$nFeature)$index[1]
  inflection <- min(nCount_inflection, nFeature_inflection)
  inflection_left <- round(inflection - inflection * 0.3)
  inflection_right <- round(inflection + inflection * 0.3)

  meta_info$RankSE <- (meta_info$nCount_rank - meta_info$nFeature_rank)^2
  if (isTRUE(fill_RankMSE)) {
    x0 <- meta_info$RankSE
    x0_nCount <- meta_info$nCount
    x0_cell <- rownames(meta_info)
    x1 <- rep(x0[-length(x0)], -diff(x0_nCount))
    x1_nCount <- rep(x0_nCount[-length(x0_nCount)], -diff(x0_nCount))
    x1_cell <- rep(x0_cell[-length(x0_cell)], -diff(x0_nCount))
    x1_cell <- paste0("Sim-", seq_len(length(x1_cell)), "@Raw-", x1_cell)
    df <- rbind(
      data.frame(RankMSE = x0, nCount = x0_nCount, cell = x0_cell),
      data.frame(RankMSE = x1, nCount = x1_nCount, cell = x1_cell)
    )
    df <- df[order(df[, "nCount"], decreasing = TRUE), ]

    for (t in seq_len(smooth_num)) {
      df[, "RankMSE"] <- runMean(df[, "RankMSE"], n = smooth_window)
      df[, "RankMSE"][is.na(df[, "RankMSE"])] <- na.omit(df[, "RankMSE"])[1]
      df[, "RankMSE"][df[, "RankMSE"] < 0] <- 0
    }
    df[, "logRankMSE"] <- log10(df[, "RankMSE"] + 1)
  } else {
    df <- meta_info
    df[, "droplets"] <- rownames(meta_info)
    df[, "RankMSE"] <- df[, "RankSE"]
    for (t in seq_len(smooth_num)) {
      df[, "RankMSE"] <- runMean(df[, "RankMSE"], n = smooth_window)
      df[, "RankMSE"][is.na(df[, "RankMSE"])] <- na.omit(df[, "RankMSE"])[1]
      df[, "RankMSE"][df[, "RankMSE"] < 0] <- 0
    }
    df[, "logRankMSE"] <- log10(df[, "RankMSE"] + 1)
  }
  # qplot(log10(1:length(df$RankMSE)), log(df$RankMSE))

  Cell_rank <- NULL
  Uncertain_rank <- NULL
  Empty_rank <- NULL
  if (isTRUE(find_rank)) {
    ## 'Cell' RankMSE valley
    df_inflection <- head(which(df[, "nCount"] <= meta_info$nCount[inflection]), 1)
    df_inflection_left <- head(which(df[, "nCount"] <= meta_info$nCount[inflection_left]), 1)
    df_inflection_right <- tail(which(df[, "nCount"] >= meta_info$nCount[inflection_right]), 1)
    crk <- df_inflection_left + find_peaks(-df$RankMSE[df_inflection_left:df_inflection_right],
      left_shoulder = df_inflection * 0.05,
      right_shoulder = df_inflection_right
    ) - 1
    # qplot(log10(1:length(df$RankMSE)), log(df$RankMSE))+geom_vline(xintercept = log10(crk))+ xlim(1, log10(df_inflection_right*1.2))
    # qplot((1:length(df$RankMSE)), log(df$RankMSE)) + xlim(1, df_inflection_right*1.2) + geom_vline(xintercept = c(crk))

    iter <- 1
    cell_max <- quantile(df[1:(crk[1] - 1), "logRankMSE"], 0.99)
    while (length(crk) > 1) {
      message("... ", length(crk), " 'Cell' valleys found. Perform automatic selection(iter=", iter, ")...")
      if (all(df[crk, "logRankMSE"] < 0)) {
        crk <- crk[length(crk)]
      }
      if (length(crk) > 1) {
        pks_diff_MSE <- diff(df[crk, "logRankMSE"])
        pks_max_MSE <- sapply(1:(length(crk) - 1), function(i) {
          quantile(df[crk[i]:crk[i + 1], "logRankMSE"], 0.99) - df[crk[i + 1], "logRankMSE"]
        })
        j <- which(pks_diff_MSE <= pks_max_MSE * tolerance | pks_max_MSE >= 0.5 * cell_max) + 1
        if (length(j) == 1) {
          j <- c(1, j)
          pks_diff_MSE <- diff(df[crk[j], "logRankMSE"])
          pks_max_MSE <- quantile(df[crk[1]:crk[2], "logRankMSE"], 0.99) - df[crk[2], "logRankMSE"]
          j <- which(pks_diff_MSE <= pks_max_MSE * tolerance | pks_max_MSE >= 0.5 * cell_max) + 1
        }
        if (length(j) == 0) {
          j <- 1
        }
        crk <- crk[j]
        iter <- iter + 1
      }
      if (length(crk) == 1) {
        message("... Automatic selection finished.")
      } else {
        message("... ", length(crk), " valleys left. Go to the next loop.")
      }
    }
    cell_count <- df[crk, "nCount"]
    Cell_rank <- max(meta_info$nCount_rank[meta_info$nCount >= cell_count])
    # qplot(log10(1:length(df$RankMSE)), log10(df$RankMSE))+geom_vline(xintercept=log10(crk))

    ## 'Empty' RankMSE valley
    maxrk <- max(which(df$nCount >= Empty_min_nCount))
    minrk <- min(crk * 5, maxrk - 10000)
    erk <- minrk + find_peaks(-df[(minrk + 1):maxrk, "RankMSE"], left_shoulder = maxrk - minrk, right_shoulder = 10000)
    erk <- erk[erk != minrk + 1]
    iter <- 1
    while (length(erk) > 1) {
      message("... ", length(erk), " 'Empty' valleys found. Perform automatic selection(iter=", iter, ")...")
      pks_diff_MSE <- diff(df[erk, "logRankMSE"])
      pks_max_MSE <- sapply(1:(length(erk) - 1), function(i) {
        quantile(df[erk[i]:erk[i + 1], "logRankMSE"], 0.99) - df[erk[i + 1], "logRankMSE"]
      })
      j <- which(-pks_diff_MSE >= pks_max_MSE * tolerance) + 1
      if (length(j) == 1) {
        j <- c(1, j)
        pks_diff_MSE <- diff(df[erk[j], "logRankMSE"])
        pks_max_MSE <- quantile(df[erk[1]:erk[2], "logRankMSE"], 0.99) - df[erk[2], "logRankMSE"]
        j <- which(-pks_diff_MSE >= pks_max_MSE * tolerance) + 1
      }
      if (length(j) == 0) {
        j <- 1
      }
      erk <- erk[j]
      iter <- iter + 1

      if (length(erk) == 1) {
        message("... Automatic selection finished.")
      } else {
        message("... ", length(erk), " valleys left. Go to the next loop.")
      }
    }
    # erk <- min(max(minrk + which.min(df[(minrk + 1):maxrk, "RankMSE"]), crk * 20), maxrk)
    empty_count <- df[erk, "nCount"]
    Empty_rank <- max(meta_info$nCount_rank[meta_info$nCount >= empty_count])
    # qplot(log10(1:length(df$RankMSE)), log10(df$RankMSE))+geom_vline(xintercept=log10(erk))

    ## 'Uncertain' RankMSE peak
    urk <- crk * 2 + which.max(df[(crk * 2 + 1):erk, "RankMSE"])
    uncertain_count <- df[urk, "nCount"]
    Uncertain_rank <- max(meta_info$nCount_rank[meta_info$nCount >= uncertain_count])
    # qplot(log10(1:length(df$RankMSE)), log10(df$RankMSE))+geom_vline(xintercept=log10(urk))
  }

  raw_df <- unique(df[, c("droplets", "RankMSE")])
  rownames(raw_df) <- raw_df$droplets
  meta_info$RankMSE <- raw_df[rownames(meta_info), "RankMSE"]
  # qplot(log10(1:length(meta_info$RankMSE)), log10(meta_info$RankMSE))+geom_vline(xintercept = log10(c(Cell_rank,Empty_rank,Uncertain_rank)))

  return(list(
    meta_info = meta_info,
    inflection = inflection,
    Cell_rank = Cell_rank,
    Uncertain_rank = Uncertain_rank,
    Empty_rank = Empty_rank
  ))
}

CalldropSplit <- function(counts, ...) {
  result <- dropSplit(counts, ...)
  return(result)
}

