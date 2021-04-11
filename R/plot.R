utils::globalVariables(c(".x", "Exp", "Feature", "Gain", "RankMSE", "nCount", "nCount_rank", "nFeature", "nFeature_rank"))

#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_histogram position_identity scale_fill_brewer labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.nCountPlot <- function(meta_info, colorBy, cell_stat_by) {
  meta_info <- subset(meta_info, nCount >= 10)
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
  p <- ggplot(meta_info, aes(x = nCount)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "nCount Histogram",
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nCount", y = "Count"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 0.5,
    )
  if (is.numeric(meta_info[, colorBy])) {
    colorBy <- "dropSplitClass"
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_histogram(
        aes(color = meta_info[, colorBy], fill = meta_info[, colorBy]),
        bins = 50, alpha = 0.5, position = position_identity()
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        scale_fill_manual(
          name = colorBy,
          values = color
        )
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_histogram(
        aes(color = meta_info[, colorBy], fill = meta_info[, colorBy]),
        bins = 50, alpha = 0.5, position = position_identity()
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        scale_fill_brewer(
          name = colorBy,
          palette = "Set1"
        )
    }
  }
  preDefinedClass <- unique(meta_info$preDefinedClass)
  p <- p + geom_vline(
    xintercept = c(
      sapply(preDefinedClass, function(x) {
        min(meta_info$nCount[meta_info$preDefinedClass == x])
      })
    ),
    color = color[preDefinedClass],
    linetype = 2
  )
}
nCountPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass") {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .nCountPlot(meta_info, colorBy, cell_stat_by)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .nCountPlot(x, colorBy, cell_stat_by)))
    }
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_histogram position_identity scale_fill_brewer labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.nFeaturePlot <- function(meta_info, colorBy, cell_stat_by) {
  meta_info <- subset(meta_info, nFeature >= 10)
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
  p <- ggplot(meta_info, aes(x = nFeature)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "nFeature Histogram",
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nFeature", y = "Count"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 0.5,
    )
  if (is.numeric(meta_info[, colorBy])) {
    colorBy <- "dropSplitClass"
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_histogram(
        aes(color = meta_info[, colorBy], fill = meta_info[, colorBy]),
        bins = 50, alpha = 0.5, position = position_identity()
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        scale_fill_manual(
          name = colorBy,
          values = color
        )
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_histogram(
        aes(color = meta_info[, colorBy], fill = meta_info[, colorBy]),
        bins = 50, alpha = 0.5, position = position_identity()
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        scale_fill_brewer(
          name = colorBy,
          palette = "Set1"
        )
    }
  }
}
nFeaturePlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass") {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .nFeaturePlot(meta_info, colorBy, cell_stat_by)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .nFeaturePlot(x, colorBy, cell_stat_by)))
    }
  }
  return(p)
}


#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_y_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.RankPlot <- function(meta_info, colorBy, cell_stat_by) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
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
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nCount_rank", y = "nFeature_rank"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (is.numeric(meta_info[, colorBy])) {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      size = 1, alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  preDefinedClass <- unique(meta_info$preDefinedClass)
  p <- p + geom_vline(
    xintercept = c(
      sapply(preDefinedClass, function(x) {
        max(meta_info$nCount_rank[meta_info$preDefinedClass == x])
      })
    ),
    color = color[preDefinedClass],
    linetype = 2
  )
  return(p)
}
RankPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass") {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .RankPlot(meta_info, colorBy, cell_stat_by)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .RankPlot(x, colorBy, cell_stat_by)))
    }
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_y_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
#' @importFrom TTR runMean
.RankMSEPlot <- function(meta_info, colorBy, cell_stat_by, ...) {
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")

  out <- RankMSE(meta_info = meta_info, ...)
  meta_info <- out$meta_info

  p <- ggplot(meta_info, aes(x = nCount_rank, y = RankMSE)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(
      title = "Smoothed RankMSE",
      subtitle = paste("#Cell:", sum(meta_info[, cell_stat_by] == "Cell")), x = "nCount_rank", y = "RankMSE"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (is.numeric(meta_info[, colorBy])) {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      size = 1,  alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = setNames(
            c("red3", "forestgreen", "steelblue", "grey80"),
            c("Cell", "Uncertain", "Empty", "Discarded")
          )
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
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
RankMSEPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass", ...) {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .RankMSEPlot(meta_info, colorBy, cell_stat_by, ...)
  if ("preDefinedClass" %in% colnames(meta_info)) {
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
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .RankMSEPlot(x, colorBy, cell_stat_by)))
    }
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.CellEntropyPlot <- function(meta_info, colorBy, cell_stat_by, highlight) {
  meta_info <- subset(meta_info, !is.na(CellEntropy))
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  cell <- intersect(highlight, rownames(meta_info))
  if (length(cell) > 0) {
    meta_info[, "highlight"] <- FALSE
    meta_info[cell, "highlight"] <- TRUE
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
  p <- ggplot(meta_info, aes(x = nCount, y = CellEntropy)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "Cell Entropy",
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nCount", y = "CellEntropy"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (is.numeric(meta_info[, colorBy])) {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      size = 1, alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(cell) > 0) {
    p <- p + geom_point(data = subset(meta_info[meta_info$highlight == TRUE, ]), color = "red", size = 1, alpha = 0.5, shape = 21)
  }
  preDefinedClass <- unique(meta_info$preDefinedClass)
  p <- p + geom_vline(
    xintercept = c(
      sapply(preDefinedClass, function(x) {
        min(meta_info$nCount[meta_info$preDefinedClass == x])
      })
    ),
    color = color[preDefinedClass],
    linetype = 2
  )
}
CellEntropyPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass", highlight = NULL) {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .CellEntropyPlot(meta_info, colorBy, cell_stat_by, highlight)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .CellEntropyPlot(x, colorBy, cell_stat_by, highlight)))
    }
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.CellRedundancyPlot <- function(meta_info, colorBy, cell_stat_by, highlight) {
  meta_info <- subset(meta_info, !is.na(CellRedundancy))
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  cell <- intersect(highlight, rownames(meta_info))
  if (length(cell) > 0) {
    meta_info[, "highlight"] <- FALSE
    meta_info[cell, "highlight"] <- TRUE
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
  p <- ggplot(meta_info, aes(x = nCount, y = CellRedundancy)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "Cell Redundancy",
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nCount", y = "CellRedundancy"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (is.numeric(meta_info[, colorBy])) {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      size = 1, alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(cell) > 0) {
    p <- p + geom_point(data = subset(meta_info[meta_info$highlight == TRUE, ]), color = "red", size = 1, alpha = 0.5, shape = 21)
  }
  preDefinedClass <- unique(meta_info$preDefinedClass)
  p <- p + geom_vline(
    xintercept = c(
      sapply(preDefinedClass, function(x) {
        min(meta_info$nCount[meta_info$preDefinedClass == x])
      })
    ),
    color = color[preDefinedClass],
    linetype = 2
  )
}
CellRedundancyPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass", highlight = NULL) {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .CellRedundancyPlot(meta_info, colorBy, cell_stat_by, highlight)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .CellRedundancyPlot(x, colorBy, cell_stat_by, highlight)))
    }
  }
  return(p)
}

#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
.CellGiniPlot <- function(meta_info, colorBy, cell_stat_by, highlight) {
  meta_info <- subset(meta_info, !is.na(CellGini))
  if (nrow(meta_info) == 0) {
    return(NULL)
  }
  cell <- intersect(highlight, rownames(meta_info))
  if (length(cell) > 0) {
    meta_info[, "highlight"] <- FALSE
    meta_info[cell, "highlight"] <- TRUE
  }
  color <- c("red3", "forestgreen", "steelblue", "grey80")
  names(color) <- c("Cell", "Uncertain", "Empty", "Discarded")
  p <- ggplot(meta_info, aes(x = nCount, y = CellGini)) +
    scale_x_log10(
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = "Cell Gini index",
      subtitle = paste0("#Cell: ", sum(meta_info[, cell_stat_by] == "Cell")),
      x = "nCount", y = "CellGini"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
  if (is.numeric(meta_info[, colorBy])) {
    p <- p + geom_point(
      aes(color = meta_info[, colorBy]),
      size = 1, alpha = 0.5, shape = 16
    ) +
      scale_color_viridis_c(
        name = colorBy
      )
  }
  if (!is.numeric(meta_info[, colorBy])) {
    if (colorBy %in% c("preDefinedClass", "dropSplitClass")) {
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_manual(
          name = colorBy,
          values = color
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      meta_info[, colorBy] <- factor(meta_info[, colorBy], levels = unique(meta_info[, colorBy]))
      p <- p + geom_point(
        aes(color = meta_info[, colorBy]),
        size = 1, alpha = 0.5, shape = 16
      ) +
        scale_color_brewer(
          name = colorBy,
          palette = "Set1"
        ) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  if (length(cell) > 0) {
    p <- p + geom_point(data = subset(meta_info[meta_info$highlight == TRUE, ]), color = "red", size = 1, alpha = 0.5, shape = 21)
  }
  preDefinedClass <- unique(meta_info$preDefinedClass)
  p <- p + geom_vline(
    xintercept = c(
      sapply(preDefinedClass, function(x) {
        min(meta_info$nCount[meta_info$preDefinedClass == x])
      })
    ),
    color = color[preDefinedClass],
    linetype = 2
  )
}
CellGiniPlot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass", highlight = NULL) {
  meta_info <- as.data.frame(meta_info)
  if ("preDefinedClass" %in% colnames(meta_info)) {
    meta_info[, "preDefinedClass"] <- factor(meta_info[, "preDefinedClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  if ("dropSplitClass" %in% colnames(meta_info)) {
    meta_info[, "dropSplitClass"] <- factor(meta_info[, "dropSplitClass"],
      levels = c("Cell", "Uncertain", "Empty", "Discarded")
    )
  }
  p <- .CellGiniPlot(meta_info, colorBy, cell_stat_by, highlight)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      meta_sp <- split.data.frame(meta_info, meta_info[[splitBy]])
      p <- c(list(Merge = p), lapply(meta_sp, function(x) .CellGiniPlot(x, colorBy, cell_stat_by, highlight)))
    }
  }
  return(p)
}

#' QC plot for droplets using meta_info generated by dropSplit.
#'
#' @param meta_info \code{importance_matrix} generated by dropSplit.
#' @param colorBy One of the column names of the \code{meta_info} or NULL. Can be a type of continuous variable or discrete variable.
#' @param splitBy One of the column names of the \code{meta_info} or NULL. Must be a type of discrete variable.
#' @param cell_stat_by One of the discrete variable column names in \code{meta_info}. Used to count the number of 'Cell' in the column.
#' @param metrics QC metrics used to plot. Possible options are:
#' \itemize{
#' \item \code{nCount}
#' \item \code{nFeature}
#' \item \code{Rank}
#' \item \code{RankMSE}
#' \item \code{CellEntropy}
#' \item \code{CellRedundancy}
#' \item \code{CellGini}
#' }
#' @param ... Other arguments passed to \code{\link{RankMSE}}.
#'
#' @return A ggplot object or a list of ggplot objects.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- dropSplit(counts)
#' pl <- QCplot(result$meta_info, metrics = "CellEntropy")
#' pl$CellEntropy$Merge
#' pl$CellEntropy$Cell
#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs guides guide_legend scale_x_log10 scale_color_brewer scale_color_manual scale_color_viridis_c annotation_logticks theme_classic theme
#' @importFrom scales trans_format math_format
#' @importFrom stats setNames
#' @importFrom methods hasArg
#' @export
QCplot <- function(meta_info, colorBy = "dropSplitClass", splitBy = "dropSplitClass", cell_stat_by = "dropSplitClass",
                   metrics = c("nCount", "nFeature", "Rank", "RankMSE", "CellEntropy", "CellRedundancy", "CellGini"), ...) {
  if (any(!metrics %in% c("nCount", "nFeature", "Rank", "RankMSE", "CellEntropy", "CellRedundancy", "CellGini"))) {
    stop("'metrics' must be in the options: 'nCount','nFeature','Rank','RankMSE','CellEntropy','CellRedundancy','CellGini'")
  }
  if (!hasArg(meta_info)) {
    stop("Parameter 'meta_info' not found.")
  }
  pl <- list()
  meta_info <- as.data.frame(meta_info)
  if (!is.null(splitBy)) {
    if (splitBy %in% colnames(meta_info)) {
      if (class(meta_info[, splitBy]) == "numeric") {
        stop("'splitBy' must be a type of discrete variable.")
      }
    } else {
      stop("'splitBy' must be one of the column names of the meta_info.")
    }
  }
  if (!is.null(colorBy)) {
    if (!colorBy %in% colnames(meta_info)) {
      stop("'colorBy' must be one of the column names of the meta_info.")
    }
  }
  args1 <- mget(names(formals()))
  args2 <- as.list(match.call())
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args1[["metrics"]] <- args1[["..."]] <- NULL

  for (metric in metrics) {
    tryCatch(expr = {
      pl[[metric]] <- suppressWarnings(base::do.call(
        what = paste0(metric, "Plot"),
        args = args1
      ))
    }, error = function(e) {
      message(e)
    })
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
#' counts <- simSimpleCounts()
#' result <- dropSplit(counts)
#' pl <- ImportancePlot(result$meta_info, result$train, result$importance_matrix, top_n = 20)
#' pl$Importance
#' pl$preDefinedClassExp
#' @importFrom ggplot2 ggplot aes geom_col geom_boxplot scale_fill_manual labs theme_classic theme
#' @importFrom stats setNames
#' @importFrom methods hasArg
#' @export
ImportancePlot <- function(meta_info, train, importance_matrix, top_n = 20) {
  if (!hasArg(meta_info)) {
    stop("Parameter 'meta_info' not found.")
  }
  if (!hasArg(train)) {
    stop("Parameter 'train' not found.")
  }
  if (!hasArg(importance_matrix)) {
    stop("Parameter 'importance_matrix' not found.")
  }
  meta_info <- as.data.frame(meta_info)
  color <- c("red3", "forestgreen", "steelblue")
  names(color) <- c("Cell", "Uncertain", "Empty")
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
      values = color
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
