#' Simulate simple counts matrix for the empty and cell-containing droplets.
#'
#' @inheritParams simComplexCounts
#'
#' @return A sparse Matrix of class "dgCMatrix".
#'
#' @examples
#' counts <- simSimpleCounts()
#' counts
#' @importFrom Matrix colSums
#' @importFrom stats rexp rgamma rpois runif
#' @export
simSimpleCounts <- function(total_gene = 30000,
                            nempty = 20000, nlarge = 2000, nsmall = 200,
                            empty_prof = NULL, empty_ngene_rate = 0.05, empty_mean_count = 100,
                            large_prof = NULL, large_ngene_rate = 0.5, large_shape = 5, large_scale = 1000,
                            small_prof = NULL, small_ngene_rate = 0.3, small_shape = 3, small_scale = 500,
                            remove_zero_drop = TRUE, remove_zero_feature = TRUE, seed = 0) {
  set.seed(seed)

  large_counts <- simCell(
    total_gene = total_gene, ncell = nlarge,
    cell_prof = large_prof, cell_ngene_rate = large_ngene_rate,
    cell_shape = large_shape, cell_scale = large_scale,
    frag = FALSE
  )
  small_counts <- simCell(
    total_gene = total_gene, ncell = nsmall,
    cell_prof = small_prof, cell_ngene_rate = small_ngene_rate,
    cell_shape = small_shape, cell_scale = small_scale,
    frag = FALSE
  )
  cell_counts <- cbind(large_counts, small_counts)
  if (is.null(empty_prof)) {
    empty_prof <- Matrix::rowSums(cell_counts)
    empty_prof <- empty_prof / sum(empty_prof)
  }
  empty_counts <- simEmpty(
    total_gene = total_gene, nempty = nempty, empty_mean_count = empty_mean_count,
    empty_prof = empty_prof, empty_ngene_rate = empty_ngene_rate
  )

  out <- cbind(large_counts, small_counts, empty_counts)
  rownames(out) <- paste0("Gene-", seq_len(total_gene))
  colnames(out) <- c(
    paste0("Cell-", "large-", "cluster1-", seq_len(ncol(large_counts))),
    paste0("Cell-", "small-", "cluster1-", seq_len(ncol(small_counts))),
    paste0("Empty-", "empty-", "cluster1-", seq_len(ncol(empty_counts)))
  )

  if (isTRUE(remove_zero_drop)) {
    out <- out[, Matrix::colSums(out) > 0]
  }
  if (isTRUE(remove_zero_feature)) {
    out <- out[Matrix::rowSums(out) > 0, ]
  }
  return(out)
}

#' Simulate a complex counts matrix including different types of empty and cell-containing droplets.
#' @description Simulation of a complex single-cell sequencing dataset.
#'
#' @param total_gene Total gene number for all the simulated counts.
#' @param disturbance A numeric value used as a weight of standard deviations when sample different distribution parameters for each cell/empty type from defined global parameters. Default is 0.2.
#' @param nempty,nlarge,nsmall Empty, large cell and small cell droplets number. If \code{remove_zero_drop} is \code{TRUE}, the final number may be smaller beacuse droplets that have all zero-valued counts will be removed.
#' @param empty_type,large_type,small_type Total number of types for \code{nempty,nlarge,nsmall}.
#' @param empty_prof,large_prof,small_prof The overall gene expression profile distribution. If provided, must be the same length with \code{total_gene}.
#' @param empty_ngene_rate,large_ngene_rate,small_ngene_rate Rate of total genes expressed in each type of droplets.
#' @param empty_mean_count Mean counts for empty droplets.
#' @param large_shape,small_shape shape parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale} control the expected mean value of counts in large or small cells.
#' @param large_scale,small_scale scale parameters in the \code{\link[stats]{GammaDist}}. \code{shape*scale^2} control the expected variance of counts in large or small cells.
#' @param large_frag,small_frag Whether simulate cell fragments from the large or small cells. Default is TRUE.
#' @param large_frag_gene,small_frag_gene Indices of cell fragment gene in profile. Default is 1:100.
#' @param large_frag_prop,small_frag_prop Proportion of the cell fragment gene counts. Default is 0.5.
#' @param remove_zero_drop Whether to remove all zero-valued droplets.
#' @param remove_zero_feature Whether to remove all zero-valued features.
#' @param seed Random seed used in simulation. Default is 0.
#'
#' @return A sparse Matrix of class "dgCMatrix".
#'
#' @examples
#' counts <- simComplexCounts()
#' counts
#' @importFrom Matrix colSums
#' @importFrom stats rexp rgamma rpois runif rnorm
#' @export
simComplexCounts <- function(total_gene = 30000, disturbance = 0.2,
                             nempty = 20000, nlarge = 5000, nsmall = 500,
                             empty_type = 5, large_type = 10, small_type = 2,
                             empty_prof = NULL, empty_ngene_rate = 0.05, empty_mean_count = 100,
                             large_prof = NULL, large_ngene_rate = 0.5, large_shape = 5, large_scale = 1000,
                             small_prof = NULL, small_ngene_rate = 0.4, small_shape = 4, small_scale = 500,
                             large_frag = TRUE, large_frag_gene = 1:50, large_frag_prop = 0.5,
                             small_frag = TRUE, small_frag_gene = 1:50, small_frag_prop = 0.5,
                             remove_zero_drop = TRUE, remove_zero_feature = TRUE, seed = 0) {
  set.seed(seed)

  large_shapes <- rnorm(n = large_type, mean = large_shape, sd = disturbance * large_shape)
  large_shapes[large_shapes < 1] <- 1
  large_scales <- rnorm(n = large_type, mean = large_scale, sd = disturbance * large_scale)
  large_scales[large_scales < 1] <- 1
  large_ngene_rates <- rnorm(n = large_type, mean = large_ngene_rate, sd = disturbance * large_ngene_rate)
  large_ngene_rates[large_ngene_rates < 0] <- 0.1
  large_ngene_rates[large_ngene_rates > 1] <- 1
  nlarges <- rnorm(n = large_type, mean = 1 / large_type, sd = disturbance / large_type)
  nlarges <- round(nlarges / sum(nlarges) * nlarge, digits = 0)
  nlarges <- nlarges[-1]
  nlarges <- c(nlarges, nlarge - sum(nlarges))
  large_frag_props <- rnorm(n = large_type, mean = large_frag_prop, sd = disturbance * large_frag_prop)
  large_frag_props[large_frag_props < 0.1] <- 0.1
  large_frag_props[large_frag_props > 0.9] <- 0.9
  large_counts_list <- list()
  large_cell_list <- list()
  large_cluster_list <- list()
  for (i in seq_len(large_type)) {
    large_counts_list[[i]] <- simCell(
      total_gene = total_gene, ncell = nlarges[i],
      cell_prof = large_prof, cell_ngene_rate = large_ngene_rates[i],
      cell_shape = large_shapes[i], cell_scale = large_scales[i],
      frag = large_frag, frag_gene = large_frag_gene, frag_gene_prop = large_frag_props[i]
    )
    large_cell_list[[i]] <- gsub(x = colnames(large_counts_list[[i]]), pattern = "-.*", replacement = "")
    large_cluster_list[[i]] <- rep(paste0("cluster", i), ncol(large_counts_list[[i]]))
  }
  large_counts <- Reduce(function(x, y) cbind(x, y), large_counts_list)
  large_cell <- Reduce(function(x, y) c(x, y), large_cell_list)
  large_cluster <- Reduce(function(x, y) c(x, y), large_cluster_list)

  small_shapes <- rnorm(n = small_type, mean = small_shape, sd = disturbance * small_shape)
  small_shapes[small_shapes < 1] <- 1
  small_scales <- rnorm(n = small_type, mean = small_scale, sd = disturbance * small_scale)
  small_scales[small_scales < 1] <- 1
  small_ngene_rates <- rnorm(n = small_type, mean = small_ngene_rate, sd = disturbance * small_ngene_rate)
  small_ngene_rates[small_ngene_rates < 0] <- 0.1
  small_ngene_rates[small_ngene_rates > 1] <- 1
  nsmalls <- rnorm(n = small_type, mean = 1 / small_type, sd = disturbance / small_type)
  nsmalls <- round(nsmalls / sum(nsmalls) * nsmall, digits = 0)
  nsmalls <- nsmalls[-1]
  nsmalls <- c(nsmalls, nsmall - sum(nsmalls))
  small_frag_props <- rnorm(n = small_type, mean = small_frag_prop, sd = disturbance * small_frag_prop)
  small_frag_props[small_frag_props < 0.1] <- 0.1
  small_frag_props[small_frag_props > 0.9] <- 0.9
  small_counts_list <- list()
  small_cell_list <- list()
  small_cluster_list <- list()
  for (i in seq_len(small_type)) {
    small_counts_list[[i]] <- simCell(
      total_gene = total_gene, ncell = nsmalls[i],
      cell_prof = small_prof, cell_ngene_rate = small_ngene_rates[i],
      cell_shape = small_shapes[i], cell_scale = small_scales[i],
      frag = small_frag, frag_gene = small_frag_gene, frag_gene_prop = small_frag_props[i]
    )
    small_cell_list[[i]] <- gsub(x = colnames(small_counts_list[[i]]), pattern = "-.*", replacement = "")
    small_cluster_list[[i]] <- rep(paste0("cluster", i), ncol(small_counts_list[[i]]))
  }
  small_counts <- Reduce(function(x, y) cbind(x, y), small_counts_list)
  small_cell <- Reduce(function(x, y) c(x, y), small_cell_list)
  small_cluster <- Reduce(function(x, y) c(x, y), small_cluster_list)

  cell_counts <- cbind(large_counts, small_counts)
  if (is.null(empty_prof)) {
    empty_prof <- Matrix::rowSums(cell_counts)
    empty_prof <- empty_prof / sum(empty_prof)
  }
  empty_mean_counts <- rnorm(n = empty_type, mean = empty_mean_count, sd = disturbance * empty_mean_count)
  empty_mean_counts[empty_mean_counts < 10] <- 10
  empty_ngene_rates <- rnorm(n = empty_type, mean = empty_ngene_rate, sd = disturbance * empty_ngene_rate)
  empty_ngene_rates[empty_ngene_rates < 0] <- 0.1
  empty_ngene_rates[empty_ngene_rates > 1] <- 1
  nemptys <- rnorm(n = empty_type, mean = 1 / empty_type, sd = 0.2 / empty_type)
  nemptys <- round(nemptys / sum(nemptys) * nempty, digits = 0)
  nemptys <- nemptys[-1]
  nemptys <- c(nemptys, nempty - sum(nemptys))
  empty_counts_list <- list()
  empty_cluster_list <- list()
  for (i in seq_len(empty_type)) {
    empty_counts_list[[i]] <- simEmpty(
      total_gene = total_gene, nempty = nemptys[i],
      empty_ngene_rate = empty_ngene_rates[i],
      empty_prof = empty_prof, empty_mean_count = empty_mean_counts[i]
    )
    empty_cluster_list[[i]] <- rep(paste0("cluster", i), ncol(empty_counts_list[[i]]))
  }
  empty_counts <- Reduce(function(x, y) cbind(x, y), empty_counts_list)
  empty_cluster <- Reduce(function(x, y) c(x, y), empty_cluster_list)

  out <- cbind(cell_counts, empty_counts)
  rownames(out) <- paste0("Gene-", seq_len(total_gene))
  colnames(out) <- c(
    paste0(large_cell, "-", "large-", large_cluster, "-", seq_len(length(large_cluster))),
    paste0(small_cell, "-", "small-", small_cluster, "-", seq_len(length(small_cluster))),
    paste0("Empty", "-", "empty-", empty_cluster, "-", seq_len(length(empty_cluster)))
  )

  if (isTRUE(remove_zero_drop)) {
    out <- out[, Matrix::colSums(out) > 0]
  }
  if (isTRUE(remove_zero_feature)) {
    out <- out[Matrix::rowSums(out) > 0, ]
  }
  return(out)
}

simEmpty <- function(total_gene = 10000, nempty = 20000, empty_mean_count = 50,
                     empty_prof = NULL, empty_ngene_rate = 0.01) {
  if (is.null(empty_prof)) {
    n <- round(total_gene * empty_ngene_rate, digits = 0)
    empty_prof <- c(
      runif(n = n, min = 0, max = 1),
      rep(0, total_gene - n, digits = 0)
    )
    empty_prof <- sample(empty_prof)
  } else {
    empty_prof <- empty_prof * sample(c(0, 1), size = length(empty_prof), replace = TRUE, prob = c(1 - empty_ngene_rate, empty_ngene_rate))
  }
  empty_prof <- empty_prof / sum(empty_prof)

  total_count <- rexp(nempty, rate = 1 / empty_mean_count)
  empty_counts <- matrix(rpois(total_gene * nempty, lambda = outer(
    empty_prof,
    total_count
  )), ncol = nempty, nrow = total_gene)
  empty_counts <- as(empty_counts, "dgCMatrix")
  colnames(empty_counts) <- paste0("Empty-", seq_len(ncol(empty_counts)))
  rownames(empty_counts) <- paste0("Gene-", seq_len(nrow(empty_counts)))
  return(empty_counts)
}

simCell <- function(total_gene = 10000, ncell = 5000,
                    cell_prof = NULL, cell_ngene_rate = 0.5,
                    cell_shape = 50, cell_scale = 50,
                    frag = TRUE, nfrag = ncell * 0.05, frag_gene = NULL,
                    frag_gene_prop = 0.8) {
  if (is.null(cell_prof)) {
    n <- ceiling(total_gene * cell_ngene_rate)
    cell_prof <- c(
      exp(rnorm(n = n, mean = 0, sd = 2)),
      rep(0, total_gene - n)
    )
    cell_prof <- sample(cell_prof)
  } else {
    cell_prof <- cell_prof * sample(c(0, 1), size = length(cell_prof), replace = TRUE, prob = c(1 - cell_ngene_rate, cell_ngene_rate))
  }
  cell_prof <- cell_prof / sum(cell_prof)

  total_count <- rgamma(ncell, shape = cell_shape, scale = cell_scale)
  cell_counts <- matrix(rpois(total_gene * ncell, lambda = outer(
    cell_prof,
    total_count
  )), ncol = ncell, nrow = total_gene)
  cell_counts <- as(cell_counts, "dgCMatrix")
  colnames(cell_counts) <- paste0("Cell-", seq_len(ncol(cell_counts)))
  rownames(cell_counts) <- paste0("Gene-", seq_len(nrow(cell_counts)))

  if (isTRUE(frag)) {
    exp_gene <- which(cell_prof > 0)
    if (is.null(frag_gene)) {
      frag_gene <- sample(exp_gene, size = 20, replace = FALSE)
    }
    nfrag_each_prop <- ceiling(nfrag / 20)
    raw_prop <- sum(cell_prof[frag_gene])
    if (raw_prop > 0 & raw_prop < frag_gene_prop) {
      for (i in seq(raw_prop + (frag_gene_prop - raw_prop) / 2, frag_gene_prop, length.out = 10)) {
        w1 <- frag_gene_prop / raw_prop
        w2 <- (1 - frag_gene_prop) / (1 - raw_prop)
        frag_prof <- cell_prof
        frag_prof[frag_gene] <- frag_prof[frag_gene] * w1
        frag_prof[-frag_gene] <- frag_prof[-frag_gene] * w2

        frag_total_count <- rexp(nfrag_each_prop, rate = 1 / (0.2 * cell_shape * cell_scale))
        frag_counts <- matrix(rpois(total_gene * nfrag_each_prop, lambda = outer(
          frag_prof,
          frag_total_count
        )), ncol = nfrag_each_prop, nrow = total_gene)
        frag_counts <- as(frag_counts, "dgCMatrix")
        colnames(frag_counts) <- paste0("Fragment-", seq_len(ncol(frag_counts)))
        rownames(frag_counts) <- paste0("Gene-", seq_len(nrow(frag_counts)))
        cell_counts <- cbind(cell_counts, frag_counts)
      }
    } else {
      message("Note: 'frag_gene' specified are zero-expressed in cell gene profile. No cell fragment generated.")
    }
  }
  return(cell_counts)
}
