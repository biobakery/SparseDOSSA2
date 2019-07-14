SparseDOSSA2 <- function(n_sample,
                         n_feature,
                         template = NULL,
                         original_feature = FALSE,
                         feature_feature_association = TRUE,
                         spike_metadata = NULL,
                         spike_strength = NULL,
                         spike_feature = NULL,
                         spike_discrete_bidirectional = TRUE,
                         biomass = NULL,
                         read_depth = 50000,
                         seed,
                         verbose = TRUE) {

  if(any(apply(feature_abd == 0, 1, all)) |
     any(apply(feature_abd == 0, 2, all))){
    warning("Feature table has all zero sample/feature (will be removed)!")
  }

  d <- nrow(feature_abd)
  n <- ncol(feature_abd)

  # read depth
  n_i <- apply(feature_abd, 2, sum)
  # Estimate MLE for pi, mu, sigma
  l_fitEM <- lapply(1:d, function(i_feature) {
    EM_theta(n_ij = feature_abd[i_feature, ],
             n_i = n_i)
  })

  # estimate F
  fit_EMFeatureParams <- Reduce("rbind",
                           lapply(l_fitEM, function(x) x$theta))
  fit_EMFeatureParams[, 3] <- log(fit_EMFeatureParams[, 3]) - log(1 - fit_EMFeatureParams[, 3])
  fit_F <- estimate_F(fit_EMFeatureParams)

  # estimate read count
  fit_readCount <- estimate_readCount(n_i)

  # generate per-feature param
  if(original_features)
    feature_paramsNew <-  fit_EMFeatureParams
  else
    feature_paramsNew <- generate_featureParam(fit_F, fit_EMFeatureParams, n_feature)

  # generate basis
  mat_basis <- generate_basis(feature_paramsNew, n_sample)

  # generate read depth
  readDepth_new <- generate_readDepth(fit_readCount, n_sample, read_depth)

  # generate read counts
  mat_count <- generate_readCount(mat_basis, readDepth_new)

  return(mat_count)
}
