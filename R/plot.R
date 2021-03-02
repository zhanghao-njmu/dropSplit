#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_y_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.RankPlot <- function(meta_info, colorBy) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  meta_info[, "nCount_rank"] <- rank(-(meta_info[, "nCount"]))
  meta_info[, "nFeature_rank"] <- rank(-(meta_info[, "nFeature"]))

  p <- ggplot(meta_info, aes(x = nCount_rank, y = nFeature_rank)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(
      title = "nCount/nFeature Rank",
      subtitle = paste0("#Cell: ", sum(meta_info$dropSplitClass == "Cell")),
      x = "nCount_rank", y = "nFeature_rank"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (class(meta_info[, colorBy]) == "numeric") {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (class(meta_info[, colorBy]) != "numeric") {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = setNames(
            c("red3", "forestgreen", "steelblue", "grey85"),
            c("Cell", "Uncertain", "Empty", "Discarded")
          )
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(unique(meta_info$preDefinedClass)) > 1) {
    p <- p + geom_vline(
      xintercept = c(
        max(meta_info$nCount_rank[meta_info$preDefinedClass == "Cell"]),
        max(meta_info$nCount_rank[meta_info$preDefinedClass == "Uncertain"]),
        max(meta_info$nCount_rank[meta_info$preDefinedClass == "Empty"])
      ),
      color = c("red3", "forestgreen", "steelblue"),
      linetype = 2
    )
  }
  return(p)
}
RankPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass") {
  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
                                           levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
                                          levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  p <- .RankPlot(meta_info, colorBy)
  if (splitBy %in% colnames(meta_info)) {
    meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
    p <- c(list(Merge = p), lapply(meta_sp, function(x) .RankPlot(x, colorBy)))
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_y_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
#' @importFrom TTR runMean
.RankMSEPlot <- function(meta_info, colorBy) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  meta_info[, "nCount_rank"] <- rank(-(meta_info[, "nCount"]))
  meta_info[, "nFeature_rank"] <- rank(-(meta_info[, "nFeature"]))
  meta_info <- meta_info[order(meta_info[, "nCount_rank"], decreasing = FALSE), ]
  meta_info[, "RankSE"] <- (meta_info[, "nCount_rank"] - meta_info[, "nFeature_rank"])^2
  x <- meta_info[, "RankSE"]
  for (t in seq_len(3)) {
    x <- runMean(x, n = 100)
    x[is.na(x)] <- na.omit(x)[1]
  }
  meta_info[, "RankMSE"] <- x
  p <- ggplot(meta_info, aes(x = nCount_rank, y = RankMSE)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(
      title = "Mean Squared Error of nCount/nFeature Rank within a window.",
      subtitle = paste("#Cell:", sum(meta_info$dropSplitClass == "Cell")), x = "nCount_rank", y = "RankMSE"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (class(meta_info[, colorBy]) == "numeric") {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (class(meta_info[, colorBy]) != "numeric") {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = setNames(
            c("red3", "forestgreen", "steelblue", "grey85"),
            c("Cell", "Uncertain", "Empty", "Discarded")
          )
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  return(p)
}
RankMSEPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass") {
  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
                                           levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
                                          levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  p <- .RankMSEPlot(meta_info, colorBy)
  p <- p + geom_vline(
    xintercept = c(
      max(meta_info$nCount_rank[meta_info$preDefinedClass == "Cell"]),
      max(meta_info$nCount_rank[meta_info$preDefinedClass == "Uncertain"]),
      max(meta_info$nCount_rank[meta_info$preDefinedClass == "Empty"])
    ),
    color = c("red3", "forestgreen", "steelblue"),
    linetype = 2
  )
  if (splitBy %in% colnames(meta_info)) {
    meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
    p <- c(list(Merge = p), lapply(meta_sp, function(x) .RankMSEPlot(x, colorBy)))
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.CellEntropyPlot <- function(meta_info, colorBy) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  meta_info <- subset(meta_info, !is.na(CellEntropy))
  p <- ggplot(meta_info, aes(x = nCount, y = CellEntropy)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "Cell Entropy",
      subtitle = paste0("#Cell: ", sum(meta_info$dropSplitClass == "Cell")),
      x = "nCount", y = "CellEntropy"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (class(meta_info[, colorBy]) == "numeric") {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (class(meta_info[, colorBy]) != "numeric") {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = setNames(
            c("red3", "forestgreen", "steelblue", "grey85"),
            c("Cell", "Uncertain", "Empty", "Discarded")
          )
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(unique(meta_info$preDefinedClass)) > 1) {
    p <- p + geom_vline(
      xintercept = c(
        min(meta_info$nCount[meta_info$preDefinedClass == "Cell"]),
        min(meta_info$nCount[meta_info$preDefinedClass == "Uncertain"]),
        min(meta_info$nCount[meta_info$preDefinedClass == "Empty"])
      ),
      color = c("red3", "forestgreen", "steelblue"),
      linetype = 2
    )
  }
}
CellEntropyPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass") {
  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
                                           levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
                                          levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  p <- .CellEntropyPlot(meta_info, colorBy)
  if (splitBy %in% colnames(meta_info)) {
    meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
    p <- c(list(Merge = p), lapply(meta_sp, function(x) .CellEntropyPlot(x, colorBy)))
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.CellGiniPlot <- function(meta_info, colorBy) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  meta_info <- subset(meta_info, !is.na(CellGini))
  p <- ggplot(meta_info, aes(x = nCount, y = CellGini)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "Cell Gini index",
      subtitle = paste0("#Cell: ", sum(meta_info$dropSplitClass == "Cell")),
      x = "nCount", y = "CellGini"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (class(meta_info[, colorBy]) == "numeric") {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (class(meta_info[, colorBy]) != "numeric") {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = setNames(
            c("red3", "forestgreen", "steelblue", "grey85"),
            c("Cell", "Uncertain", "Empty", "Discarded")
          )
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(unique(meta_info$preDefinedClass)) > 1) {
    p <- p + geom_vline(
      xintercept = c(
        min(meta_info$nCount[meta_info$preDefinedClass == "Cell"]),
        min(meta_info$nCount[meta_info$preDefinedClass == "Uncertain"]),
        min(meta_info$nCount[meta_info$preDefinedClass == "Empty"])
      ),
      color = c("red3", "forestgreen", "steelblue"),
      linetype = 2
    )
  }
}
CellGiniPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass") {
  meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
                                           levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
                                          levels = c("Cell", "Uncertain", "Empty", "Discarded")
  )
  p <- .CellGiniPlot(meta_info, colorBy)
  if (splitBy %in% colnames(meta_info)) {
    meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
    p <- c(list(Merge = p), lapply(meta_sp, function(x) .CellGiniPlot(x, colorBy)))
  }
  return(p)
}

#' QC plot for droplets using meta_info generated by dropSplit.
#'
#' @param meta_info \code{importance_matrix} generated by dropSplit.
#' @param colorBy One of the column names of the \code{meta_info}. Can be a type of continuous variable or discrete variable.
#' @param splitBy One of the column names of the \code{meta_info}. Must be a type of discrete variable.
#' @param metrics QC metrics used to plot. Possible options are:
#' \itemize{
#' \item \code{Rank}
#' \item \code{RankMSE}
#' \item \code{CellEntropy}
#' \item \code{CellGini}
#' }
#'
#' @return A ggplot object or a list of ggplot objects.
#'
#' @examples
#' counts <- simCounts()
#' result <- dropSplit(counts)
#' pl <- QCPlot(result$meta_info, metrics = "CellEntropy")
#' pl$CellEntropy$Merge
#' pl$CellEntropy$Cell
#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
#' @export
QCPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass",
                   metrics = c("Rank", "RankMSE", "CellEntropy", "CellGini")) {
  if (!any(metrics %in% c("Rank", "RankMSE", "CellEntropy", "CellGini"))) {
    stop("'metrics' must be at least one of the 'Rank','RankMSE','CellEntropy','CellGini' ")
  }
  pl <- list()
  for (metric in metrics) {
    if (!metric %in% c("Rank", "RankMSE", "CellEntropy", "CellGini")) {
      warning("metric:", metric, " is not valid. Skipped.", immediate. = TRUE)
    } else {
      args <- list(meta_info = meta_info, colorBy = colorBy, splitBy = splitBy)

      tryCatch(expr = {
        pl[[metric]] <- suppressWarnings(base::do.call(
          what = paste0(metric, "Plot"),
          args = args
        ))
      }, error = function(e) {
        message(e)
      })
    }
  }
  return(pl)
}

#' Plot the feature importance
#'
#' @param meta_info \code{importance_matrix} generated by dropSplit.
#' @param train \code{train} dataset generated by dropSplit.
#' @param importance_matrix \code{importance_matrix} generated by dropSplit.
#' @param top_n Number of the top important features to plot.
#'
#' @return A a list of ggplot objects.
#'
#' @examples
#' counts <- simCounts()
#' colnames(counts) <- paste0("Cell-", 1:ncol(counts))
#' result <- dropSplit(counts)
#' pl <- ImportancePlot(result$meta_info, result$train, result$importance_matrix, top_n = 20)
#' pl$Importance
#' pl$preDefinedClassExp
#' @importFrom ggplot2 ggplot aes geom_col geom_boxplot scale_fill_manual labs theme_classic theme
#' @importFrom stats setNames
#' @export
ImportancePlot <- function(meta_info, train, importance_matrix, top_n = 20) {
  df <- importance_matrix[1:min(top_n, nrow(importance_matrix)), ]
  df <- df[order(df$Gain, decreasing = FALSE), ]
  df$Feature <- factor(df$Feature, levels = df$Feature)
  p1 <- ggplot(df, aes(y = Feature, x = Gain)) +
    geom_col(color = "black", fill = "steelblue") +
    labs(
      title = "Feature Importance",
      x = "Gain", y = "Feature"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  barcodes <- rownames(meta_info)
  preDefinedClass <- meta_info$preDefinedClass
  dropSplitClass <- meta_info$dropSplitClass

  index <- barcodes %in% rownames(train)
  barcodes <- barcodes[index]
  preDefinedClass <- preDefinedClass[index]
  dropSplitClass <- dropSplitClass[index]
  train <- train[barcodes, df$Feature]
  dat <- data.frame(
    Exp = as.numeric(train),
    Feature = rep(as.character(df$Feature), each = nrow(train)),
    preDefinedClass = rep(as.character(preDefinedClass), times = nrow(df)),
    dropSplitClass = rep(as.character(dropSplitClass), times = nrow(df))
  )
  dat$Feature <- factor(dat$Feature, levels = df$Feature)
  dat$preDefinedClass <- factor(dat$preDefinedClass, levels = c("Cell", "Empty"))
  dat$dropSplitClass <- factor(dat$dropSplitClass, levels = c("Cell", "Uncertain", "Empty", "Discarded"))
  p2 <- ggplot(dat, aes(y = Feature, x = log2(Exp), fill = preDefinedClass)) +
    geom_boxplot(outlier.size = 0) +
    scale_fill_manual(
      name = "preDefinedClass",
      values = setNames(
        c("red3", "steelblue"),
        c("Cell", "Empty")
      )
    ) +
    labs(
      title = "Feature Expression in droplets group by preDefinedClass",
      x = "log2(Expression)", y = "Feature"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  p3 <- ggplot(dat, aes(y = Feature, x = log2(Exp), fill = dropSplitClass)) +
    geom_boxplot(outlier.size = 0) +
    scale_fill_manual(
      name = "dropSplitClass",
      values = setNames(
        c("red3", "forestgreen", "steelblue"),
        c("Cell", "Uncertain", "Empty")
      )
    ) +
    labs(
      title = "Feature Expression in droplets group by dropSplitClass",
      x = "log2(Expression)", y = "Feature"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  return(list(
    Importance = p1,
    preDefinedClassExp = p2,
    dropSplitClassExp = p3
  ))
}
