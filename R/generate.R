generate_featureParam <- function(fit_F, feature_params, n_feature) {
  if(fit_F$pi0_sigma2 > 0)
    sigma2_ind <- rbinom(n = n_feature, size = 1, prob = 1 - fit_F$pi0_sigma2)
  else 
    sigma2_ind <- rep(1, length = n_feature)
  
  # simulate parameters for non-zero sigma2s
  feature_paramsNonZeroSigma2 <- feature_params[feature_params[, 2] > 0, ]
  nonzeroSigma2_ind <- sigma2_ind[sigma2_ind > 0]
  mat_simulatedParams <- matrix(NA, nrow = sum(nonzeroSigma2_ind),
                                ncol = 3)
  while(TRUE) {
    nFeature_nonzero <- sum(nonzeroSigma2_ind)
    sample_mixture <- feature_paramsNonZeroSigma2[sample.int(n = nrow(feature_paramsNonZeroSigma2),
                                                             size = nFeature_nonzero,
                                                             replace = TRUE), ]
    sample_difference <- mvtnorm::rmvnorm(n = nFeature_nonzero, sigma = fit_F$K_nonzero)
    samples_nonzero <- sample_mixture + sample_difference
    mat_simulatedParams[nonzeroSigma2_ind == 1] <- samples_nonzero
    nonzeroSigma2_ind <- (mat_simulatedParams[, 2] <= 0) * 1
    if(sum(nonzeroSigma2_ind) == 0) break
  }
  
  nFeature_zero <- sum(sigma2_ind == 0)
  if(nFeature_zero > 0) {
    feature_paramsZeroSigma2 <- feature_params[feature_params[, 2] == 0, ]
    sample_mixture <- feature_paramsZeroSigma2[sample.int(n = nrow(feature_paramsZeroSigma2),
                                                             size = nFeature_zero,
                                                             replace = TRUE), ]
    sample_difference <-  mvtnorm::rmvnorm(n = nFeature_zero, sigma = fit_F$K_zero)
    sample_mixture[, c(1, 3)] <- sample_mixture[, c(1, 3)] + sample_difference
    mat_simulatedParams <- rbind(mat_simulatedParams, sample_mixture)
  }
  return(mat_simulatedParams)
}

generate_basis <- function(feature_params, n_sample) {
  feature_params[, 3] <- exp(feature_params[, 3]) / (1 + exp(feature_params[, 3]))
  mat_basis <- sapply(1:nrow(feature_params), function(i_feature) {
    generate_ZILogitNM(feature_params[i_feature, 1], 
                       feature_params[i_feature, 2],
                       feature_params[i_feature, 3],
                       n_sample = n_sample)
  })
  mat_basis <- apply(mat_basis, 1, function(x) x / sum(x))
}

generate_ZILogitNM <- function(mu, sigma2, pi0, n_sample) {
  samples_nonzero <- rnorm(n = n_sample, mean = mu, sd = sqrt(sigma2))
  samples_nonzero <- exp(samples_nonzero) / (1 + exp(samples_nonzero))
  samples_ind <- rbinom(n = n_sample, size = 1, p = 1 - pi0)
  return(samples_nonzero * samples_ind)
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
