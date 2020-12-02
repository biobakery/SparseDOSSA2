generate_a <- function(n, 
                       feature_param, Omega,
                       maxit = 10, verbose = FALSE) {
  i_iter <- 0
  samples_a <- matrix(NA, 
                      nrow = nrow(feature_param),
                      ncol = 0)
  rownames(samples_a) <- rownames(feature_param)
  
  while(TRUE) {
    if(i_iter + 1 > maxit)
      stop("Can't satisfy conditions!")
    i_iter <- i_iter + 1
    if(verbose) 
      print(i_iter)
    
    samples_a <- cbind(samples_a,
                       t(rcopulasso(n = n,
                                    pi0 = feature_param[, "pi0"],
                                    mu = feature_param[, "mu"],
                                    sigma = feature_param[, "sigma"],
                                    Omega = Omega)))
    
    ind_nonzero <- apply(samples_a > 0, 2, any) ## FIXME??
    if(sum(ind_nonzero) >= n) {
      samples_a <- samples_a[, ind_nonzero][, seq_len(n)]
      colnames(samples_a) <- paste0("Sample", seq_len(n))
      return(samples_a)
    }
  }
}

rcopulasso <- function(n, pi0, mu, sigma, Omega) {
  if(length(pi0) != length(mu) |
     length(pi0) != length(sigma) |
     length(pi0) != nrow(Omega))
    stop("Parameter dimensions must agree!")
  
  # sample marginals
  mat_amarginals <- 
    vapply(seq_len(length(pi0)),
           function(i)
             rZILogN_one(n = n,
                         pi0 = pi0[i],
                         mu = mu[i],
                         sigma = sigma[i]),
           rep(0.0, n))
  # arrange from smallest to largest for shuffling
  mat_amarginals <- 
    apply(mat_amarginals, 2, function(x) x[order(x)])
  
  # sample ranks
  mat_rank <- 
    mvtnorm::rmvnorm(n = n, sigma = solve(Omega))
  mat_rank <- apply(mat_rank, 2, rank)
  
  mat_a <- 
    vapply(seq_len(length(pi0)),
           function(i)
             mat_amarginals[, i, drop = TRUE][mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  
  return(mat_a)
}

rZILogN_one <- function(n, pi0, mu, sigma) {
    return(exp(rnorm(n = n, mean = mu, sd = sigma)) *
             rbinom(n = n, size = 1, prob = 1 - pi0))
}

generate_featureParam <- function(F_fit, param_original, n_feature) {
  if(!all(names(F_fit) == c("p0_sigma", "K_nonzero", "K_zero")))
    stop("F_fit is not of the correct format!")
  
  # transform pi parameter
  param_original[, "pi0"] <- logit(param_original[, "pi0"])
  
  if(F_fit$p0_sigma > 0)
    ind_sigma <- rbinom(n = n_feature, size = 1, 
                        prob = 1 - F_fit$p0_sigma) == 1
  else
    ind_sigma <- rep(TRUE, length = n_feature)
  
  # simulate parameters for non-zero sigmas
  param_original_nonZeroSigma <-
    param_original[param_original[, "sigma"] > 0, , drop = FALSE]
  ind_sigma_nonzero_tosimulate <- ind_sigma[ind_sigma]
  param_simulated <- 
    matrix(NA, nrow = sum(ind_sigma_nonzero_tosimulate), ncol = 3)
  dimnames(param_simulated) <- list(paste0("Feature", seq_len(n_feature)),
                                    c("pi0", "mu", "sigma"))
  while(TRUE) {
    nFeature_nonzero <- sum(ind_sigma_nonzero_tosimulate)
    simulated_mixture <- 
      param_original_nonZeroSigma[sample.int(
        n = nrow(param_original_nonZeroSigma),
        size = nFeature_nonzero,
        replace = TRUE), , drop = FALSE]
    simulated_difference <- mvtnorm::rmvnorm(n = nFeature_nonzero, 
                                             sigma = F_fit$K_nonzero)
    simulated_mixture <- simulated_mixture + simulated_difference
    param_simulated[ind_sigma_nonzero_tosimulate, ] <- simulated_mixture
    ind_sigma_nonzero_tosimulate <- param_simulated[, "sigma"] <= 0
    if(!any(ind_sigma_nonzero_tosimulate)) break
  }
  
  nFeature_zero <- sum(!ind_sigma)
  if(nFeature_zero > 0) {
    param_original_zeroSigma <- 
      param_original[param_original[, "sigma"] == 0, , drop = FALSE]
    simulated_mixture <- 
      param_original_zeroSigma[sample.int(
        n = nrow(param_original_zeroSigma),
        size = nFeature_zero,
        replace = TRUE), , drop = FALSE]
    simulated_difference <-  mvtnorm::rmvnorm(n = nFeature_zero, 
                                              sigma = F_fit$K_zero)
    simulated_mixture[, c("pi0", "mu")] <- 
      simulated_mixture[, c("pi0", "mu"), drop = FALSE] +
      simulated_difference
    param_simulated <- rbind(param_simulated, simulated_mixture)
  }
  
  param_simulated[, "pi0"] <- expit(param_simulated[, "pi0"])
  return(param_simulated)
}


#' @importFrom magrittr %>%
generate_feature_metadata_spike_df <- 
  function(features,
           perc_feature_spiked_metadata,
           n_metadata,
           effect_size,
           spike_metadata = c("both", 
                              "abundance",
                              "prevalence")) {
    n_feature_spiked <- ceiling(length(features) *
                                  perc_feature_spiked_metadata)
    if(spike_metadata == "both")
      spike_metadata <- c("abundance", "prevalence")
    
    feature_metadata_spike_df <- 
      purrr::map2_dfr(seq_len(n_metadata),
                      effect_size,
                      function(metadatum_i, effect_size_i) {
                        data.frame(metadata_datum = metadatum_i,
                                   effect_size = 
                                     effect_size_i * 
                                     sample(c(-1, 1),
                                            size = n_feature_spiked,
                                            replace = TRUE),
                                   feature_spiked = 
                                     sample(features,
                                            size = n_feature_spiked,
                                            replace = FALSE),
                                   stringsAsFactors = FALSE)
                      })
    feature_metadata_spike_df <- 
      tidyr::expand_grid(feature_metadata_spike_df,
                         associated_property = spike_metadata) %>% 
      dplyr::select(metadata_datum, 
                    feature_spiked,
                    associated_property,
                    effect_size) %>% 
      dplyr::arrange(metadata_datum,
                     associated_property,
                     -effect_size) %>% 
      as.data.frame()
    
    return(feature_metadata_spike_df)
  }

generate_depth <- function(mu_depth, 
                           sigma_depth, 
                           n, median_depth) {
  depth <- exp(rnorm(n = n, mean = mu_depth, sd = sigma_depth))
  depth <- round(depth / median(depth) * median_depth)
  return(depth)
}

generate_count <- function(rel, depth) {
  mat_count <- 
    vapply(seq_along(depth), 
           function(i_sample) 
             rmultinom(n = 1, size = depth[i_sample], prob = rel[, i_sample]),
           rep(0, nrow(rel)))
  dimnames(mat_count) <- dimnames(rel)
  return(mat_count)
}  
