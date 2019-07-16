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
    ind_sigma2 <- rbinom(n = n_feature, size = 1, prob = 1 - fit_F$p0_sigma2)
  else
    ind_sigma2 <- rep(1, length = n_feature)

  # simulate parameters for non-zero sigma2s
  feature_paramNonzeroSigma2 <- feature_param[feature_param[, 2] > 0, ]
  ind_sigma2_nonzero <- ind_sigma2[ind_sigma2 > 0]
  mat_simulatedParam <- matrix(NA, nrow = sum(ind_sigma2_nonzero),
                               ncol = 3)
  while(TRUE) {
    nFeature_nonzero <- sum(ind_sigma2_nonzero)
    simulated_mixture <- feature_paramNonzeroSigma2[sample.int(n = nrow(feature_paramNonzeroSigma2),
                                                               size = nFeature_nonzero,
                                                               replace = TRUE), ]
    simulated_difference <- mvtnorm::rmvnorm(n = nFeature_nonzero, sigma = fit_F$K_nonzero)
    mat_simulatedNonzero <- simulated_mixture + simulated_difference
    mat_simulatedParam[ind_sigma2_nonzero == 1] <- mat_simulatedNonzero
    ind_sigma2_nonzero <- (mat_simulatedParam[, 2] <= 0) * 1
    if(sum(ind_sigma2_nonzero) == 0) break
  }

  nFeature_zero <- sum(ind_sigma2 == 0)
  if(nFeature_zero > 0) {
    feature_paramZeroSigma2 <- feature_param[feature_param[, 2] == 0, ]
    simulated_mixture <- feature_paramZeroSigma2[sample.int(n = nrow(feature_paramZeroSigma2),
                                                            size = nFeature_zero,
                                                            replace = TRUE), ]
    simulated_difference <-  mvtnorm::rmvnorm(n = nFeature_zero, sigma = fit_F$K_zero)
    simulated_mixture[, c(1, 3)] <- simulated_mixture[, c(1, 3)] + simulated_difference
    mat_simulatedParam <- rbind(mat_simulatedParam, simulated_mixture)
  }

  mat_simulatedParam[, 3] <- exp(mat_simulatedParam[, 3]) / (1 + exp(mat_simulatedParam[, 3]))
  return(mat_simulatedParam)
}

generate_basis <- function(feature_param, n_sample, fit_C, seed) {
  set.seed(seed)
  mat_basis_original <- sapply(1:nrow(feature_param), function(i_feature) {
    generate_ZILogitNM(feature_param[i_feature, 1],
                       feature_param[i_feature, 2],
                       feature_param[i_feature, 3],
                       n_sample = n_sample)
  })

  if(!is.null(fit_C)) {
    simulated_rank <- copula::rCopula(fit_C, n = nrow(feature_param))
    simulated_rank <- simulated_rank %>%
      apply(2, rank)
    mat_basis_shuffled <- sapply(1:ncol(sample),
                                 function(j) sort(sample[, j])[sample_rank[, j]])
  } else {
    mat_basis_shuffled <- mat_basis_original
  }

  mat_basis <- apply(mat_basis, 1, function(x) x / sum(x))
}

generate_ZILogitN <- function(mu, sigma2, pi, n_sample) {
  simulated_nonzero <- rnorm(n = n_sample, mean = mu, sd = sqrt(sigma2))
  simulated_nonzero <- exp(simulated_nonzero) / (1 + exp(simulated_nonzero))
  ind_nonzero <- rbinom(n = n_sample, size = 1, p = pi)
  return(simulated_nonzero * ind_nonzero)
}

generate_spiked_basis <- function(feature_param, n_sample) {
  feature_param[, 3] <- exp(feature_param[, 3]) / (1 + exp(feature_param[, 3]))
  mat_basis <- sapply(1:nrow(feature_param), function(i_feature) {
    generate_ZILogitNM(feature_param[i_feature, 1],
                       feature_param[i_feature, 2],
                       feature_param[i_feature, 3],
                       n_sample = n_sample)
  })
  mat_basis <- apply(mat_basis, 1, function(x) x / sum(x))
}

generate_ZILogitN <- function(mu, sigma2, pi0, n_sample) {
  mat_simulatedNonzero <- rnorm(n = n_sample, mean = mu, sd = sqrt(sigma2))
  mat_simulatedNonzero <- exp(mat_simulatedNonzero) / (1 + exp(mat_simulatedNonzero))
  samples_ind <- rbinom(n = n_sample, size = 1, p = 1 - pi0)
  return(mat_simulatedNonzero * samples_ind)
}

generate_readDepth <- function(LN_params, n_sample, read_depth) {
  samples <- exp(rnorm(n = n_sample, mean = LN_params[1], sd = sqrt(LN_params[2])))
  samples <- round(samples / mean(samples) * read_depth)
  return(samples)
}

generate_readCount <- function(mat_p, n_i) {
  sapply(1:length(n_i), function(i_sample) {
    rmultinom(n = 1, size = n_i[i_sample], prob = mat_p[, i_sample])
  })
}
