SparseDOSSA2 <- function(n_sample = NULL,
                         n_feature = NULL,
                         template = NULL,
                         same_features = FALSE,
                         feature_feature_association = TRUE,
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

  # set control parameters if specified

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

    # estimate F
    if(!same_features) {
      if(verbose) message("same_features is FALSE, estimating distribution of per-feature parameters...")
      fit_F <- estimate_F(dfFit_featureParam, control = control$kd)
    }

    # estimate correlation copula
    if(feature_feature_association) {
      if(verbose) message("feature_feature_association is TRUE, estimating feature-feature associations...")
      mat_posterior <- Reduce("rbind", lapply(lFit_featureParam, function(x) x$hidden_param$mu_posterior_ij))
      mat_posterior <- exp(mat_posterior) / (1 + exp(mat_posterior))
      mat_posterior[template == 0] <- 0
      fit_C <- estimate_C(mat_posterior, control = control$copula, seed = seed)
    }

    # estimate read count
    if(verbose) message("Estimating read count distribution...")
    fit_readCount <- estimate_readCount(n_i)

    if(verbose) message("Parameter estimation from template dataset complete!")
  }


    # generate per-feature param

    if(is.null(n_feature))
      featureParam_new <-  dfFit_featureParam
    else {
      featureParam_new <- generate_featureParam(fit_F, fit_EMFeatureParams, n_feature)
    }



    # estimate C


    # generate basis
    mat_basis <- generate_basis(feature_paramsNew, n_sample)


  if(!is.null(spike_metadata)) {

  }



  # fit_EMFeatureParams[, 3] <- log(fit_EMFeatureParams[, 3]) - log(1 - fit_EMFeatureParams[, 3])






  # generate read depth
  readDepth_new <- generate_readDepth(fit_readCount, n_sample, read_depth)

  # generate read counts
  mat_count <- generate_readCount(mat_basis, readDepth_new)

  return(mat_count)
}
