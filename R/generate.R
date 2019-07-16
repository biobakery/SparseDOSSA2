#' Generate per-feature parameters (mu, sigma, pi) given pre-fit non-parameter density kernels
#'
#' @param fit_F a list (fit generated fro estimate_F)
#' @param feature_param data frame of estimated per-feature parameters (in order mu, sigma2, and pi)
#' @param n_feature number of features to generate
#' @param seed random seed
#'
#' @return matrix of per-feature parameters (in order mu, sigma2, and pi)
generate_featureParam <- function(fit_F, feature_param, n_feature, seed) {
  set.seed(seed)

  if(!all(names(fit_F) == c("p0_sigma2", "K_nonzero", "K_zero")))
    stop("fit_F is not of the correct format!")

  # transform pi parameter
  feature_param[, 3] <- log(feature_param[, 3]) - log(1 - feature_param[, 3])

  if(fit_F$p0_sigma2 > 0)
    ind_sigma2 <- rbinom(n = n_feature, size = 1, prob = 1 - fit_F$p0_sigma2) == 1
  else
    ind_sigma2 <- rep(TRUE, length = n_feature)

  # simulate parameters for non-zero sigma2s
  feature_paramNonzeroSigma2 <- feature_param[feature_param[, 2] > 0, , drop = FALSE]
  ind_sigma2_nonzero_tosimulate <- ind_sigma2[ind_sigma2]
  mat_simulatedParam <- matrix(NA, nrow = sum(ind_sigma2_nonzero_tosimulate),
                               ncol = 3)
  while(TRUE) {
    nFeature_nonzero <- sum(ind_sigma2_nonzero_tosimulate)
    simulated_mixture <- feature_paramNonzeroSigma2[sample.int(n = nrow(feature_paramNonzeroSigma2),
                                                               size = nFeature_nonzero,
                                                               replace = TRUE), , drop = FALSE]
    simulated_difference <- mvtnorm::rmvnorm(n = nFeature_nonzero, sigma = fit_F$K_nonzero)
    mat_simulatedNonzero <- simulated_mixture + simulated_difference
    mat_simulatedParam[ind_sigma2_nonzero_tosimulate, ] <- mat_simulatedNonzero
    ind_sigma2_nonzero_tosimulate <- mat_simulatedParam[, 2] <= 0
    if(!any(ind_sigma2_nonzero_tosimulate)) break
  }

  nFeature_zero <- sum(!ind_sigma2)
  if(nFeature_zero > 0) {
    feature_paramZeroSigma2 <- feature_param[feature_param[, 2] == 0, , drop = FALSE]
    simulated_mixture <- feature_paramZeroSigma2[sample.int(n = nrow(feature_paramZeroSigma2),
                                                            size = nFeature_zero,
                                                            replace = TRUE), , drop = FALSE]
    simulated_difference <-  mvtnorm::rmvnorm(n = nFeature_zero, sigma = fit_F$K_zero)
    simulated_mixture[, c(1, 3)] <- simulated_mixture[, c(1, 3), drop = FALSE] +
      simulated_difference
    mat_simulatedParam <- rbind(mat_simulatedParam, simulated_mixture)
  }

  mat_simulatedParam[, 3] <- exp(mat_simulatedParam[, 3]) / (1 + exp(mat_simulatedParam[, 3]))
  return(mat_simulatedParam)
}

#' Generate basis (i.e., pij) matrix, with the option to incorporate copula shuffling for correlation structure
#'
#' @param feature_param data frame/matrix of generated (if same_features = FALSE) or estimated
#' (if same_feature = TRUE) per-feature parameters (in order mu, sigma2, and pi)
#' @param n_sample number of samples to generate
#' @param fit_C a copula fit (fit generated from estimate_C). If null then no shuffling will be performed.
#' @param seed random seed
#'
#' @return a list of matrices, in order:
#' mat_basis_original is the originally generated matrix
#' mat_basis_shuffled is the copula shuffled matrix (is equal to mat_basis_original if fit_C is null)
#' mat_basis is the normalized basis matrix for proper multinomial sampling
generate_basis <- function(feature_param, n_sample, fit_C = NULL,
                           seed) {
  set.seed(seed)

  # keep sampling until no samples have all-zero values
  ind_sample_tosimulate <- rep(TRUE, n_sample)
  tMat_basis_original <- tMat_basis_shuffled <- matrix(NA_real_, nrow = n_sample, ncol = nrow(feature_param))
  while(TRUE) {
    simulated_tMat_basis_original <- sapply(1:nrow(feature_param), function(i_feature) {
      generate_ZILogitN(feature_param[i_feature, 1],
                        feature_param[i_feature, 2],
                        feature_param[i_feature, 3],
                        n_sample = n_sample)
    })

    if(!is.null(fit_C)) {
      if(nrow(feature_param) != fit_C@dimension)
        stop("Number of features disagree in feature_param and fit_C!")
      tSimulated_rank <- copula::rCopula(fit_C, n = n_sample)
      tSimulated_rank <- apply(tSimulated_rank, 2, rank)
      simulated_tMat_basis_shuffled <- sapply(1:ncol(simulated_tMat_basis_original),
                                              function(i_feature)
                                                sort(simulated_tMat_basis_original[, i_feature])[
                                                  tSimulated_rank[, i_feature]])
    } else {
      simulated_tMat_basis_shuffled <- simulated_tMat_basis_original
    }

    ind_simulated_nonzero <- apply(simulated_tMat_basis_shuffled > 0, 1, any)
    n_simulated <- min(sum(ind_sample_tosimulate), sum(ind_simulated_nonzero))
    tMat_basis_original[ind_sample_tosimulate, ][1:n_simulated, ] <-
      simulated_tMat_basis_original[ind_simulated_nonzero, , drop = FALSE][1:n_simulated, ]
    tMat_basis_shuffled[ind_sample_tosimulate, ][1:n_simulated, ] <-
      simulated_tMat_basis_shuffled[ind_simulated_nonzero, , drop = FALSE][1:n_simulated, ]
    ind_sample_tosimulate[ind_sample_tosimulate][1:n_simulated] <- FALSE

    if(!any(ind_sample_tosimulate)) break
  }

  mat_basis <- apply(tMat_basis_shuffled, 1, function(x) x / sum(x))
  colnames(tMat_basis_original) <-
    colnames(tMat_basis_shuffled) <-
    rownames(mat_basis) <- rownames(feature_param)
  rownames(tMat_basis_original) <-
    rownames(tMat_basis_shuffled) <-
    colnames(mat_basis) <- paste0("Sample", 1:n_sample)

  return(list(mat_basis_original = t(tMat_basis_original),
              mat_basis_shuffled = t(tMat_basis_shuffled),
              mat_basis = mat_basis))
}

#' Generate zero-inflated logistic normal per-feature basis (pij for fixed j)
#'
#' @param mu mean parameter for the normal component
#' @param sigma2 variance parameter for the normal component
#' @param pi prevalence (i.e. 1 - zero proportion) parameter
#' @param n_sample number of samples to generate
#'
#' @return vector of zero-infated proportions for a given feature
generate_ZILogitN <- function(mu, sigma2, pi, n_sample) {
  simulated_nonzero <- rnorm(n = n_sample, mean = mu, sd = sqrt(sigma2))
  simulated_nonzero <- exp(simulated_nonzero) / (1 + exp(simulated_nonzero))
  ind_nonzero <- rbinom(n = n_sample, size = 1, p = pi)
  return(simulated_nonzero * ind_nonzero)
}

# generate_spiked_basis <- function(feature_param, n_sample) {
#   feature_param[, 3] <- exp(feature_param[, 3]) / (1 + exp(feature_param[, 3]))
#   mat_basis <- sapply(1:nrow(feature_param), function(i_feature) {
#     generate_ZILogitN(feature_param[i_feature, 1],
#                       feature_param[i_feature, 2],
#                       feature_param[i_feature, 3],
#                       n_sample = n_sample)
#   })
#   mat_basis <- apply(mat_basis, 1, function(x) x / sum(x))
# }

#' Generate log-normal per-sample read counts
#'
#' @param mu mean parameter
#' @param sigma2 variance parameter
#' @param n_sample number of samples to generate
#' @param read_depth median read depth to normalize to
#' @param seed random seed
#'
#' @return vector of total read depths
generate_readDepth <- function(mu, sigma2, n_sample, read_depth, seed) {
  set.seed(seed)
  samples <- exp(rnorm(n = n_sample, mean = mu, sd = sqrt(sigma2)))
  samples <- round(samples / median(samples) * read_depth)
  return(samples)
}

#' Generate multinomial sampled read counts matrix
#'
#' @param mat_p feature x sample matrix of underlying relative abundance
#' @param n_i vector of total read count across samples
#'
#' @return feature x sample matrix of microbial read count
generate_readCount <- function(mat_p, n_i, seed) {
  set.seed(seed)
  mat_count <- sapply(1:length(n_i), function(i_sample) {
    rmultinom(n = 1, size = n_i[i_sample], prob = mat_p[, i_sample])
  })
  dimnames(mat_count) <- dimnames(mat_p)
  return(mat_count)
}
