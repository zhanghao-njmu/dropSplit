# CellCalling function -----------------------------------------------

#' Call cell-containg droplets using different methods.
#'
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.
#' @param method A cell-calling method from the following options:
#' \itemize{
#' \item \link[=dropSplit]{dropSplit}
#' \item \link[=CallCellRangerV2]{CellRangerV2}
#' \item \link[=CallCellRangerV3]{CellRangerV3}
#' \item \link[=CallEmptyDrops]{EmptyDrops}
#' \item \link[=CallzUMIs]{zUMIs}
#' }
#' @param seed Random seed used in the function. Default is 0.
#' @param ... Arguments passed to the cell-calling method.
#' @return A list including at least \code{DataFrame} named 'meta_info' which contains the final classification column named '\code{method}Class'.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- CellCalling(counts, method = "dropSplit")
#' head(result)
#' @importFrom methods hasArg
#' @export
CellCalling <- function(counts, method = "dropSplit", seed = 0, ...) {
  if (!hasArg(counts)) {
    stop("'counts' must be specified.")
  }
  if (!method %in% c("dropSplit", "CellRangerV2", "CellRangerV3", "EmptyDrops", "zUMIs")) {
    stop("'method' must be one of 'dropSplit','CellRangerV2','CellRangerV3','EmptyDrops','zUMIs'")
  }
  if (length(colnames(counts)) != ncol(counts) | length(rownames(counts)) != nrow(counts)) {
    stop("'counts' matrix must have both row(feature) names and column(cell) names.")
  }
  if (any(Matrix::colSums(counts) <= 0)) {
    stop("'counts' has droplets that nCount<=0.")
  }
  args1 <- mget(names(formals()))
  args2 <- as.list(match.call())
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args1[["..."]] <- args1[["method"]] <- args1[["seed"]] <- NULL

  set.seed(seed)

  message("\n############# Start CellCalling using ", method, " #############\n")
  tryCatch(expr = {
    out <- suppressWarnings(base::do.call(
      what = paste0("Call", method),
      args = args1
    ))
  }, error = function(e) {
    message(e)
    stop("Stop CellCalling.")
  })
  if (class(out) != "list") {
    result <- list(meta_info = out)
  } else {
    result <- out
  }
  message("\n############# CellCalling Finished #############\n")

  return(result)
}

#' EmptyDrops method used to recognize cell-containing droplets.
#'
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.
#' @param FDR Maximum FDR value to call a droplet as non-empty(labeled as 'Cell'). Default is 0.01.
#' @param lower A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets. Default is 100.
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations. Default is 10000.
#' @param ... Other arguments passed to \code{\link[DropletUtils]{emptyDrops}}.
#' @return A \code{DataFrame} containing the classification column named 'EmptyDropsClass'.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- CallEmptyDrops(counts)
#' head(result)
#' @importFrom DropletUtils emptyDrops
#' @importFrom S4Vectors DataFrame
#' @export
CallEmptyDrops <- function(counts, FDR = 0.01, lower = 100, niters = 10000, ...) {
  meta_info <- emptyDrops(counts, lower = lower, niters = niters, ...)
  EmptyDropsClass <- ifelse(meta_info$FDR <= FDR, "Cell", "Empty")
  EmptyDropsClass[is.na(EmptyDropsClass)] <- "Empty"
  nlimited_droplets <- sum(table(Sig = EmptyDropsClass, Limited = meta_info$Limited)["Empty", "TRUE"])
  if (nlimited_droplets > 0) {
    message(">>> ", nlimited_droplets, " 'Empty' droplets have a limited p-value. A lower p-value could be obtained by increasing niters.")
  }
  meta_info[, "EmptyDropsClass"] <- EmptyDropsClass
  meta_info <- DataFrame(meta_info)

  message(">>> EmptyDrops identified ", sum(meta_info$EmptyDropsClass == "Cell"), " cell-containing droplets.")
  return(meta_info)
}


#' zUMI method used to recognize cell-containing droplets.
#'
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.
#' @return A \code{DataFrame} containing the classification column named 'zUMIsClass'.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- CallzUMIs(counts)
#' head(result)
#' @importFrom inflection uik
#' @importFrom S4Vectors DataFrame
#' @export
CallzUMIs <- function(counts) {
  meta_info <- data.frame(row.names = colnames(counts))
  meta_info$nCount <- Matrix::colSums(counts)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  raw_cell_order <- rownames(meta_info)
  meta_info <- meta_info[order(meta_info$nCount_rank, decreasing = FALSE), ]
  meta_info$nCount_cumsum <- cumsum(meta_info$nCount)

  if (length(unique(as.integer(quantile(meta_info[, "nCount_rank"])))) != 5) {
    meta_info[nrow(meta_info), "nCount_rank"] <- meta_info[nrow(meta_info), "nCount_rank"] + 1
  }
  ntop <- floor(uik(x = meta_info[, "nCount_rank"], y = meta_info[, "nCount_cumsum"]))
  meta_info[, "zUMIsClass"] <- "Empty"
  meta_info[1:ntop, "zUMIsClass"] <- "Cell"
  meta_info <- meta_info[raw_cell_order, ]
  meta_info <- DataFrame(meta_info)

  message(">>> zUMIs identified ", sum(meta_info$zUMIsClass == "Cell"), " cell-containing droplets.")
  return(meta_info)
}


#' CellRanger_v2 method used to recognize cell-containing droplets.
#'
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.
#' @param recovered_cells Expected number of recovered cells. Default is 3000.
#' @param recovered_cells_quantile Quantile of the top \code{recovered_cells} barcodes by total UMI counts. Default is 0.99.
#' @return A \code{DataFrame} containing the classification column named 'CellRangerV2Class'.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- CallCellRangerV2(counts)
#' head(result)
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @export
CallCellRangerV2 <- function(counts, recovered_cells = 3000, recovered_cells_quantile = 0.99) {
  meta_info <- data.frame(row.names = colnames(counts))
  meta_info$nCount <- Matrix::colSums(counts)
  meta_info$nCount_rank <- rank(-(meta_info$nCount))
  filtered_indices <- filter_cellular_barcodes_ordmag(meta_info$nCount, recovered_cells, recovered_cells_quantile)
  meta_info[, "CellRangerV2Class"] <- "Empty"
  meta_info[filtered_indices, "CellRangerV2Class"] <- "Cell"
  meta_info <- DataFrame(meta_info)

  message(">>> CellRangerV2 identified ", sum(meta_info$CellRangerV2Class == "Cell"), " cell-containing droplets.")
  return(meta_info)
}


#' CellRanger_v3 method used to recognize cell-containing droplets.
#'
#' @param counts A \code{matrix} object or a \code{dgCMatrix} object which columns represent features and rows represent droplets.
#' @inheritParams CallCellRangerV2
#' @param n_candidate_barcodes Number of additional barcodes to consider after the initial cell calling. Default is 20000.
#' @param n_partitions Number of partitions (max number of barcodes to consider for ambient estimation). Default is 90000.
#' @param min_umis_nonambient Minimum number of UMIS per barcode to consider after the initial cell calling. Default is 500.
#' @param min_umi_frac_of_median Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling. Default is 0.01.
#' @param max_adj_pvalue Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling. Default is 0.01.
#' @return A \code{DataFrame} containing the classification column named 'CellRangerV3Class'.
#'
#' @examples
#' counts <- simSimpleCounts()
#' result <- CallCellRangerV3(counts)
#' head(result)
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @export
CallCellRangerV3 <- function(counts, recovered_cells = 3000, recovered_cells_quantile = 0.99,
                             n_candidate_barcodes = 20000, n_partitions = 90000,
                             min_umis_nonambient = 500, min_umi_frac_of_median = 0.01, max_adj_pvalue = 0.01) {
  message(">>> Use 'CallCellRangerV2' to perform initial cell calling...")
  meta_info <- CallCellRangerV2(counts, recovered_cells, recovered_cells_quantile)
  orig_cell_bcs <- rownames(meta_info)[which(meta_info$CellRangerV2Class == "Cell")]
  message(">>> Identify low RNA content cells from remaining droplets...")

  cr3 <- find_nonambient_barcodes(
    matrix = counts, orig_cell_bcs = orig_cell_bcs,
    n_partitions = n_partitions,
    n_candidate_barcodes = n_candidate_barcodes,
    min_umi_frac_of_median = min_umi_frac_of_median,
    min_umis_nonambient = min_umis_nonambient,
    max_adj_pvalue = max_adj_pvalue
  )
  meta_info[, "CellRangerV3Class"] <- meta_info[, "CellRangerV2Class"]
  if (nrow(cr3) > 0) {
    meta_info[cr3$barcode, "pvalues_adj"] <- cr3$pvalues_adj
    meta_info[cr3$barcode, "eval_bcs"] <- cr3$eval_bcs
    meta_info[cr3$barcode, "CellRangerV3Class"] <- ifelse(cr3$is_nonambient, "Cell", "Empty")
  }
  meta_info <- DataFrame(meta_info)
  message("An additional ", sum(cr3$is_nonambien), " cells were identified.")

  message("\n>>> CellRangerV3 identified ", sum(meta_info$CellRangerV3Class == "Cell"), " cell-containing droplets.")
  return(meta_info)
}


# CellRanger cell_calling.py -----------------------------------------------

#' Estimate a gene expression profile by Simple Good Turing
#' @param matrix Sparse matrix of all counts
#' @param barcode_indices Barcode indices to use
#' @param nz_feat Indices of features that are non-zero at least once
#'
#' @return Estimated probabilities of length len(nz_feat).
estimate_profile_sgt <- function(matrix, barcode_indices, nz_feat) {

  # Initial profile estimate
  prof_mat <- matrix[, barcode_indices]

  profile <- Matrix::rowSums(prof_mat[nz_feat, ])
  zero_feat <- which(profile == 0)

  # Simple Good Turing estimate
  # p_smoothed, p0
  res <- sgt_proportions(profile[which(profile != 0)])
  p_smoothed <- res[["pstar"]]
  p0 <- res[["p0"]]

  # Distribute p0 equally among the zero elements.
  p0_i <- p0 / length(zero_feat)

  profile_p <- rep(p0_i, length(nz_feat))
  profile_p[which(profile != 0)] <- p_smoothed

  if (!isTRUE(all.equal(sum(profile_p), 1))) {
    stop("sum(profile_p): ", sum(profile_p), " != 1")
  }
  return(profile_p)
}


#' Estimate a gene expression profile on a given subset of barcodes.
#' @description Use Good-Turing to smooth the estimated profile.
#' @param matrix Sparse matrix of all counts.
#' @param use_bcs Indices of barcodes to use.
#' @return A list including:
#' \itemize{
#' \item \code{use_feats} use_features
#' \item \code{bg_profile_p} Estimated probabilities of length use_features.
#' }
#'
est_background_profile_sgt <- function(matrix, use_bcs) {

  # Use features that are nonzero anywhere in the data
  use_feats <- which(Matrix::rowSums(matrix) != 0)

  # Estimate background profile
  bg_profile_p <- estimate_profile_sgt(matrix, use_bcs, use_feats)

  return(list(use_feats = use_feats, bg_profile_p = bg_profile_p))
}


#' Call barcodes as being sufficiently distinct from the ambient profile
#'
#' @param matrix Full expression matrix.
#' @param orig_cell_bcs Strings of initially-called cell barcodes.
#' @inheritParams CallCellRangerV3
#' @return A data.frame.
#' @importFrom stats p.adjust
find_nonambient_barcodes <- function(matrix, orig_cell_bcs,
                                     n_partitions = 90000,
                                     n_candidate_barcodes = 20000,
                                     min_umis_nonambient = 500,
                                     min_umi_frac_of_median = 0.01,
                                     max_adj_pvalue = 0.01) {

  # Estimate an ambient RNA profile
  umis_per_bc <- Matrix::colSums(matrix)
  bc_order <- order(umis_per_bc)

  # Take what we expect to be the barcodes associated with empty partitions.
  if (n_partitions > length(bc_order)) {
    warning("n_partitions(", n_partitions, ") is larger than the number of barcodes(", length(bc_order), "). Reset n_partitions value to ", length(bc_order), ".", immediate. = TRUE)
    n_partitions <- length(bc_order)
  }
  empty_bcs <- rev(bc_order)[((n_partitions / 2) + 1):n_partitions]
  empty_bcs <- sort(empty_bcs)

  # Require non-zero barcodes
  nz_bcs <- which(umis_per_bc != 0)
  nz_bcs <- sort(nz_bcs)

  use_bcs <- intersect(empty_bcs, nz_bcs)

  if (length(use_bcs) > 0) {
    res <- est_background_profile_sgt(matrix, use_bcs)
    eval_features <- res[["use_feats"]]
    ambient_profile_p <- res[["bg_profile_p"]]
  } else {
    eval_features <- c()
    ambient_profile_p <- c()
  }

  # Choose candidate cell barcodes
  orig_cell_bc_set <- unique(orig_cell_bcs)
  orig_cells <- colnames(matrix) %in% orig_cell_bc_set

  # Look at non-cell barcodes above a minimum UMI count
  eval_bcs <- seq_len(ncol(matrix))
  names(eval_bcs) <- colnames(matrix)
  eval_bcs[orig_cells] <- NA

  median_initial_umis <- median(umis_per_bc[orig_cells])
  min_umis <- as.integer(max(min_umis_nonambient, round(ceiling(median_initial_umis * min_umi_frac_of_median))))
  message("Median UMIs of initial cell calls: ", median_initial_umis)
  message("Min UMIs of initial cell calls: ", min(umis_per_bc[orig_cells]))
  message("Min UMIs of droplets used for prediction: ", min_umis)

  eval_bcs[umis_per_bc < min_umis] <- NA
  n_unmasked_bcs <- length(eval_bcs) - sum(is.na(eval_bcs))
  if (n_unmasked_bcs == 0) {
    warning("No barcodes left. \nYou may decrease the 'min_umis_nonambient' value if want to distinguish more cells from ambient droplets.")
  }

  umis_per_bc_mask <- umis_per_bc
  umis_per_bc_mask[which(is.na(eval_bcs))] <- NA
  # Take the top n_candidate_barcodes by UMI count, of barcodes that pass the above criteria
  # tail(umis_per_bc_mask[order(umis_per_bc_mask)[1:n_unmasked_bcs]],100)
  eval_bcs <- order(umis_per_bc_mask)[1:n_unmasked_bcs]
  names(eval_bcs) <- names(umis_per_bc_mask)[order(umis_per_bc_mask)[1:n_unmasked_bcs]]
  eval_bcs <- eval_bcs[max(length(eval_bcs) - n_candidate_barcodes + 1, 1):length(eval_bcs)]

  message("Number of candidate droplets used for prediction: ", length(eval_bcs))
  message("Range candidate droplets UMIs: ", min(umis_per_bc[eval_bcs]), ",", max(umis_per_bc[eval_bcs]))

  eval_mat <- matrix[eval_features, eval_bcs]

  # Compute observed log-likelihood of barcodes being generated from ambient RNA
  obs_loglk <- eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

  # Simulate log likelihoods
  res2 <- simulate_multinomial_loglikelihoods(
    profile_p = ambient_profile_p,
    umis_per_bc = umis_per_bc[eval_bcs],
    num_sims = 10000, verbose = TRUE
  )
  distinct_ns <- res2[["distinct_n"]]
  sim_loglk <- res2[["loglk"]]

  # Compute p-values
  pvalues <- compute_ambient_pvalues(
    umis_per_bc = umis_per_bc[eval_bcs],
    obs_loglk = obs_loglk,
    sim_n = distinct_ns,
    sim_loglk = sim_loglk
  )

  pvalues_adj <- p.adjust(pvalues, method = "BH")
  is_nonambient <- pvalues_adj <= max_adj_pvalue

  result <- data.frame(
    barcode = colnames(matrix)[eval_bcs],
    eval_bcs = eval_bcs,
    log_likelihood = obs_loglk,
    pvalues = pvalues,
    pvalues_adj = pvalues_adj,
    is_nonambient = is_nonambient
  )
  return(result)
}


# CellRanger stat.py -------------------------------------------------------
#' Compute the multinomial log PMF for many barcodes
#'
#' @param matrix Matrix of UMI counts (feature x barcode)
#' @param profile_p Multinomial probability vector
#'
#' @return Log-likelihood for each barcode
#' @importFrom stats dmultinom
eval_multinomial_loglikelihoods <- function(matrix, profile_p) {
  loglk <- rep(0, ncol(matrix))
  chunk_size <- 100000
  index <- 1:ncol(matrix)
  index_chunk <- split(index, ceiling(seq_along(index) / chunk_size))
  message("Split matrix into ", length(index_chunk), " chunks:")
  for (i in 1:length(index_chunk)) {
    message("... Calculate log PMF for chunk: ", i)
    chunk <- index_chunk[[i]]
    loglk[chunk] <- apply(matrix[, chunk], 2, function(x) {
      dmultinom(x, size = sum(x), prob = profile_p, log = TRUE)
    })
  }
  return(loglk)
}

#' Simulate draws from a multinomial distribution for various values of N.
#' @description    Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )
#'
#' @param profile_p Probability of observing each feature.
#' @param umis_per_bc UMI counts per barcode (multinomial N).
#' @param num_sims Number of simulations per distinct N value.
#' @param jump Vectorize the sampling if the gap between two distinct Ns exceeds this.
#' @param n_sample_feature_block Vectorize this many feature samplings at a time.
#' @param verbose verbose. Default is FALSE.
#' @return A list including:
#' \itemize{
#' \item \code{distinct_ns} an array containing the distinct N values that were simulated.
#' \item \code{log_likelihoods} a length(distinct_ns) x num_sims matrix containing the simulated log likelihoods.
#' }
#' @importFrom stats rmultinom dmultinom
#'
simulate_multinomial_loglikelihoods <- function(profile_p, umis_per_bc,
                                                num_sims = 1000, jump = 1000,
                                                n_sample_feature_block = 1000000, verbose = FALSE) {
  distinct_n <- sort(unique(umis_per_bc))
  loglk <- matrix(rep(0, length(distinct_n) * num_sims), ncol = num_sims)
  num_all_n <- max(distinct_n) - min(distinct_n)
  if (verbose) {
    message("Number of distinct N supplied: ", length(distinct_n))
    message("Range of N: ", num_all_n)
    message("Number of features: ", length(profile_p))
  }
  sampled_features <- sample(length(profile_p), size = n_sample_feature_block, prob = profile_p, replace = TRUE)
  k <- 1

  log_profile_p <- log(profile_p)

  for (sim_idx in seq_len(num_sims)) {
    curr_counts <- c(rmultinom(n = 1, prob = profile_p, size = distinct_n[1]))
    # eval_multinomial_loglikelihoods(matrix = curr_counts,)
    curr_loglk <- dmultinom(x = curr_counts, size = distinct_n[1], prob = profile_p, log = TRUE)
    loglk[1, sim_idx] <- curr_loglk

    for (i in seq(2, length(distinct_n))) {
      step <- distinct_n[i] - distinct_n[i - 1]
      if (step >= jump) {
        # Instead of iterating for each n, sample the intermediate ns all at once
        curr_counts <- curr_counts + rmultinom(n = 1, prob = profile_p, size = step)
        curr_loglk <- dmultinom(x = curr_counts, size = distinct_n[i], prob = profile_p, log = TRUE)
      } else {
        # Iteratively sample between the two distinct values of n
        for (n in seq(distinct_n[i - 1] + 1, distinct_n[i])) {
          j <- sampled_features[k]
          k <- k + 1
          if (k >= n_sample_feature_block) {
            # Amortize this operation
            sampled_features <- sample(length(profile_p), size = n_sample_feature_block, prob = profile_p, replace = TRUE)
            k <- 1
            curr_counts[j] <- curr_counts[j] + 1
            curr_loglk <- curr_loglk + log_profile_p[j] + log(n / curr_counts[j])
          }
        }
      }
      loglk[i, sim_idx] <- curr_loglk
    }
  }
  return(list(distinct_n = distinct_n, loglk = loglk))
}


#' Compute p-values for observed multinomial log-likelihoods
#'
#' @param umis_per_bc UMI counts per barcode
#' @param obs_loglk  Observed log-likelihoods of each barcode deriving from an ambient profile
#' @param sim_n  Multinomial N for simulated log-likelihoods
#' @param sim_loglk Simulated log-likelihoods of maxtix. nrow=length(sim_n), ncol=num_simulations
#'
#' @return pvalues
#'
compute_ambient_pvalues <- function(umis_per_bc, obs_loglk, sim_n, sim_loglk) {
  if (length(umis_per_bc) != length(obs_loglk)) {
    stop("length(umis_per_bc) != length(obs_loglk)")
  }
  if (nrow(sim_loglk) != length(sim_n)) {
    stop("nrow(sim_loglk) != length(sim_n)")
  }

  # Find the index of the simulated N for each barcode
  sim_n_idx <- sapply(umis_per_bc, function(x) {
    if (x <= min(sim_n)) {
      return(1)
    } else {
      if (x >= max(sim_n)) {
        return(length(sim_n))
      } else {
        return(max(which(x > sim_n)) + 1)
      }
    }
  })
  num_sims <- ncol(sim_loglk)

  num_barcodes <- length(umis_per_bc)

  pvalues <- rep(0, num_barcodes)

  for (i in seq_len(num_barcodes)) {
    # message(i)
    num_lower_loglk <- sum(sim_loglk[sim_n_idx[i], ] < obs_loglk[i])
    pvalues[i] <- (1 + num_lower_loglk) / (1 + num_sims)
  }
  return(pvalues)
}


find_within_ordmag <- function(x, baseline_idx) {
  x_ascending <- sort(x)
  baseline <- x_ascending[length(x_ascending) - baseline_idx + 1]
  cutoff <- max(1, round(0.1 * baseline))
  # Return the index corresponding to the cutoff in descending order
  index <- max(which(cutoff >= x_ascending)) + 1
  return(length(x) - index)
}

filter_cellular_barcodes_ordmag <- function(bc_counts, recovered_cells = 3000, recovered_cells_quantile = 0.99) {
  max_filtered_bcs <- recovered_cells * 6
  nonzero_bc_counts <- bc_counts[bc_counts > 0]
  baseline_bc_idx <- round(recovered_cells) * (1 - recovered_cells_quantile)
  baseline_bc_idx <- min(baseline_bc_idx, length(nonzero_bc_counts) - 1)

  # Bootstrap sampling; run algo with many random samples of the data
  top_n_boot <- c()
  for (i in seq_len(100)) {
    top_n_boot[i] <- find_within_ordmag(
      x = sample(x = nonzero_bc_counts, size = length(nonzero_bc_counts)),
      baseline_idx = baseline_bc_idx
    )
  }

  # Get the filtered barcodes
  top_n_bcs_mean <- mean(top_n_boot)
  top_n <- round(top_n_bcs_mean)
  top_bc_idx <- sort(order(bc_counts, decreasing = TRUE)[1:top_n])
  return(top_bc_idx = top_bc_idx)
}

# CellRanger sgt.py ---------------------------------------------------------------------

# Simple Good-Turing estimator.
# Based on S implementation in William A. Gale & Geoffrey Sampson (1995) Good-turing frequency estimation without tears,
# Journal of Quantitative Linguistics, 2:3, 217-237, DOI: 10.1080/09296179508590051

.averaging_transform <- function(r, nr) {
  d <- c(1, diff(r))
  dr <- c(
    0.5 * (d[2:length(d)] + d[1:(length(d) - 1)]),
    d[length(d)]
  )
  return(nr / dr)
}

.rstest <- function(r, coef) {
  return(r * (1 + 1 / r)^(1 + coef))
}

#' Make a Simple Good-Turing estimate of the frequencies.
#'
#' @param xr Non-zero item frequencies
#' @param xnr Non-zero frequencies of frequencies
#' @return A list including:
#' \itemize{
#' \item \code{rstar} The adjusted non-zero frequencies
#' \item \code{p0} The total probability of unobserved items
#' }
#' @importFrom stats lm
simple_good_turing <- function(xr, xnr) {
  xN <- sum(xr * xnr)
  xnrz <- .averaging_transform(xr, xnr)
  fit <- lm(log(xnrz) ~ log(xr))
  intercept <- as.numeric(fit$coefficients[[1]])
  slope <- as.numeric(fit$coefficients[[2]])
  if (slope > -1) {
    message("The log-log slope is > -1 (", slope, "); the SGT estimator is not applicable to these data.")
  }
  xrst <- .rstest(xr, slope)
  xrstrel <- xrst / xr
  xrtry <- xr == c(xr[2:length(xr)] - 1, 0)
  xrstarel <- rep(0, length(xr))
  xrstarel[xrtry] <- (xr[xrtry] + 1) / xr[xrtry] * c(xnr[2:length(xnr)], 0)[xrtry] / xnr[xrtry]


  # Determine when to switch from GT to LGT estimates
  tursd <- rep(1, length(xr))

  for (i in seq_len(length(xr))) {
    if (isTRUE(xrtry[i])) {
      tursd[i] <- (i + 2) / xnr[i] * sqrt(xnr[i + 1] * (1 + xnr[i + 1] / xnr[i]))
    }
  }
  xrstcmbrel <- rep(0, length(xr))
  useturing <- TRUE
  for (r in seq_len(length(xr))) {
    if (!useturing) {
      xrstcmbrel[r] <- xrstrel[r]
    } else {
      if (abs(c(xrstrel[r] - xrstarel[r])) * (1 + r) / tursd[r] > 1.65) {
        xrstcmbrel[r] <- xrstarel[r]
      } else {
        useturing <- FALSE
        xrstcmbrel[r] <- xrstrel[r]
      }
    }
  }

  # Renormalize the probabilities for observed objects
  sumpraw <- sum(xrstcmbrel * xr * xnr / xN)

  xrstcmbrel <- xrstcmbrel * (1 - xnr[1] / xN) / sumpraw
  p0 <- xnr[1] / xN

  return(list(rstar = xr * xrstcmbrel, p0 = p0))
}


#' Use Simple Good-Turing estimate to adjust for unobserved items
#'
#' @param frequencies Nonzero frequencies of items
#' @return A list including:
#' \itemize{
#' \item \code{pstar} The adjusted non-zero proportions
#' \item \code{p0} The total probability of unobserved items
#' }
sgt_proportions <- function(frequencies) {
  if (0 %in% frequencies) {
    stop("frequencies cannot contain zero value.")
  }
  a <- table(frequencies)
  freqfreqs <- rep(0, max(frequencies) + 1)
  names(freqfreqs) <- 0:max(frequencies)
  freqfreqs[names(a)] <- a
  use_freqs <- as.numeric(names(freqfreqs[freqfreqs != 0]))
  freqfreqs <- as.numeric(freqfreqs)

  if (length(use_freqs) < 10) {
    stop("Too few non-zero frequency items (", length(use_freqs), "). Aborting SGT.")
  }

  xr <- use_freqs
  xnr <- freqfreqs[use_freqs + 1]
  res <- simple_good_turing(xr, xnr)

  # rstar contains the smoothed frequencies.
  # Map each original frequency r to its smoothed rstar.
  rstar <- res[["rstar"]]
  p0 <- res[["p0"]]
  rstar_dict <- setNames(object = rstar, nm = use_freqs)
  rstar_sum <- sum(xnr * rstar)
  rstar_i <- rstar_dict[as.character(frequencies)]
  pstar <- (1 - p0) * (rstar_i / rstar_sum)

  if (!isTRUE(all.equal(p0 + sum(pstar), 1))) {
    stop("p0 + sum(pstar): ", (p0 + sum(pstar)), " != 1")
  }
  return(list(pstar = pstar, p0 = p0))
}


# frequencies <- use_freqs <- xr <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 36, 41, 43, 45, 46, 47, 50, 71, 84, 101, 105, 121, 124, 146, 162, 193, 199, 224, 226, 254, 257, 339, 421, 456, 481, 483, 1140, 1256, 1322, 1530, 2131, 2395, 6925, 7846)
# xnr <- c(120, 40, 24, 13, 15, 5, 11, 2, 2, 1, 3, 2, 1, 1, 3, 1, 3, 2, 3, 3, 3, 2, 2, 1, 2, 2, 1, 2, 2, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# reticulate::repl_python()
