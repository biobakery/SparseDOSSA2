#' Sparse Data Observations for Simulating Synthetic Abundance v2
#'
#' @param n_sample number of samples to simulate
#' @param n_feature number of features to simulate
#' @param template template dataset to learn community structures from
#' @param same_features should SparseDOSSA2 generate random feautures, or the same ones as
#' in the template dataset?
#' @param feature_feature_association only valid when same_feature = TRUE. Should SparseDOSSA2 additionally
#' model the feature-feature association structure in the template dataset?
#' @param spike_metadata
#' @param spike_strength
#' @param spike_feature
#' @param spike_discrete_bidirectional
#' @param biomass
#' @param read_depth target median sample read depths of the simulated microbial counts
#' @param seed random seed
#' @param control nested list of control parameters for the various model fitting components of SparseDOSSA2
#' @param verbose
#'
#' @return
#' @export
SparseDOSSA2 <- function(n_sample = NULL,
                         n_feature = NULL,
                         template = NULL,
                         same_features = FALSE,
                         feature_feature_association = FALSE,
                         spike_metadata = NULL,
                         spike_strength = NULL,
                         spike_feature = NULL,
                         spike_discrete_bidirectional = TRUE,
                         biomass = NULL,
                         read_depth = 50000,
                         seed,
                         control = list(em = list(maxiter.outer = 1000,
                                                  maxiter.inner = 1000,
                                                  reltol.outer = 1e-8,
                                                  factr.inner = 1e7),
                                        kd = list(estimator = "Hscv"),
                                        copula = list(method = "spearman",
                                                      random = TRUE,
                                                      R = 50)),
                         verbose = TRUE) {
  # parameter checking

  # set controls if specified

  set.seed(seed)

  # estimate feature parameters from template if provided
  if(!is.null(template)) {
    if(verbose) message("Template count table is provided and will be used to estimate parameters!")

    # some sanity check to ensure is count table
    template <- as.matrix(template)
    if(!all(template %% 1 == 0))
      stop("Template appears to not be count table!")

    check.features <- apply(template == 0, 1, all)
    check.samples <- apply(template == 0, 2, all)
    if(any(check.features) | any(check.samples)){
      warning("Template count table has all zero sample/feature (will be removed)!",)
      template <- template[!check.features, !check.samples, drop = FALSE]
      warning(nrow(template), " features and ", ncol(template), "samples remain.")
    }

    # estimate per-feautre parameters
    if(verbose) message("Estimating per-feature parameters...")
    ni_template <- apply(template, 2, sum)
    # Estimate MLE for pi, mu, sigma
    lFit_featureParam <- lapply(1:nrow(template), function(i_feature) {
      estimate_featureParam(y_ij = template[i_feature, ],
                            n_i = ni_template,
                            control = control$em)
    })
    dfFit_featureParam <- Reduce("rbind",
                                 lapply(lFit_featureParam, function(x) x$theta))
    rownames(dfFit_featureParam) <- rownames(template)

    # estimate F
    if(!same_features) {
      if(verbose) message("same_features is FALSE, estimating distribution of per-feature parameters...")
      fit_F <- estimate_F(feature_param = dfFit_featureParam,
                          control = control$kd)
    } else fit_F <- NULL

    # estimate correlation copula
    if(same_features & feature_feature_association) {
      if(verbose) message("Both same_features and feature_feature_association are TRUE, ",
                          "estimating feature-feature associations...")
      mat_posterior <- Reduce("rbind", lapply(lFit_featureParam, function(x) x$hidden_param$mu_posterior_ij))
      mat_posterior <- exp(mat_posterior) / (1 + exp(mat_posterior))
      mat_posterior[template == 0] <- 0
      fit_C <- estimate_C(feature_abd = mat_posterior,
                          seed = seed,
                          control = control$copula)
    } else fit_C <- NULL

    # estimate read count
    if(verbose) message("Estimating read count distribution...")
    fit_readCount <- estimate_readCount(n_i = n_i)

    if(verbose) message("Parameter estimation from template dataset complete!")
  }


  # generate per-feature params
  if(verbose) message("Generating per-feature parameters...")
  if(same_features)
    featureParam_new <-  dfFit_featureParam
  else {
    featureParam_new <- generate_featureParam(fit_F = fit_F,
                                              feature_param = dfFit_featureParam,
                                              n_feature = n_feature,
                                              seed = seed)
    rownames(featureParam_new) <- paste0("Feature", 1:n_feature)
  }

  # generate basis matrix
  if(verbose) message("Generating basis relative abundance matrix...")
  lMat_basis <- generate_basis(feature_param = featureParam_new,
                               n_sample = n_sample,
                               fit_C = fit_C,
                               seed = seed)


  # if(!is.null(spike_metadata)) {
  #
  # }

  if(verbose) message("Generating count matrix...")
  # generate read depth
  readDepth_new <- generate_readDepth(mu = fit_readCount["mu"],
                                      sigma2 = fit_readCount["sigma2"],
                                      n_sample = n_sample,
                                      read_depth = read_depth,
                                      seed = seed)
  # generate read counts
  mat_count <- generate_readCount(mat_p = lMat_basis$mat_basis,
                                  n_i = readDepth_new,
                                  seed = seed)

  # aggregate results
  results <- list(counts = mat_count,
                  template_fits = list(em = lFit_featureParam,
                                       kd = fit_F,
                                       copula = fit_C,
                                       depth = fit_readCount),
                  hidden_params = list(feature_param = featureParam_new,
                                       mats_basis = lMat_basis),
                  control = control)

  return(results)
}
