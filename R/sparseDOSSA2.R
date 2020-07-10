#' SparseDOSSA2 model fitting to microbiome datasets
#'
#' @param data n x p matrix of microbial count observations
#' @param lambdas vector of sparsity tuning parameters
#' @param K number of cross-validation folds
#' @param control list of control parameters
#'
#' @return
#' @export
#'
#' @examples
SparseDOSSA2_fit <- function(data, 
                             lambdas,
                             K = 5,
                             control = list()) {
  control <- do.call(control_EM, control)
  if(!is.null(control$debug_dir))
    dir.create(control$debug_dir)
  
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  if(control$verbose) message("Performing full data fit")
  control_tmp <- control
  control_tmp$debug_dir <- paste0(control_tmp$debug_dir, "K0/")
  l_fits_full <- future::future({
    EM(data = data,
       lambdas = lambdas,
       control = control_tmp)$l_fits
  })
  l_fits_full <- future::value(l_fits_full)
  
  # CV fits
  CV_folds <- make_CVfolds(n = nrow(data), K = K)
  if(!is.null(control$debug_dir))
    save(CV_folds, file = paste0(control$debug_dir, "CV_folds.RData"))
  l_results_CV <- list()
  for(k in seq_len(K)) {
    if(control$verbose) message("Performing CV K=", k)
    l_results_CV[[k]] <- future::future({
      data_training <- data[CV_folds != k, ]
      data_testing <- data[CV_folds == k, ]
      control_tmp <- control
      control_tmp$debug_dir <- paste0(control$debug_dir, "K", k, "/")
      
      result_k <- EM(data = data_training,
                     lambdas = lambdas,
                     control = control_tmp)
      
      # Fill in parameters estimates for features not present in training data
      for(i_lambda in seq_along(lambdas))
        result_k$l_fits[[i_lambda]]$fit <- 
        fill_estimates_CV(result_k$l_fits[[i_lambda]]$fit,
                          l_fits_full[[i_lambda]]$fit,
                          result_k$l_filtering$ind_feature)
      
      # Calculate ll in testing data
      l_logLik <- future.apply::future_lapply(
        seq_along(lambdas),
        function(i_lambda) {
          params <- result_k$l_fits[[i_lambda]]$fit
          future.apply::future_sapply(
            seq_len(nrow(data_testing)),
            function(i_sample) {
              logLik <- 
                dx(x = data_testing[i_sample, , drop = TRUE],
                   pi0 = params$pi0, mu = params$mu, sigma = params$sigma, 
                   Omega = params$Omega, Sigma = params$Sigma,
                   control = control$control_numint,
                   log.p = TRUE)
            })
        })
      
      list(l_fits = result_k$l_fits,
           l_logLik = l_logLik)
    })
  }
  l_results_CV <- future::values(l_results_CV)
  
  # Aggregate CV results
  logLik_CV <- 
    sapply(seq_along(lambdas),
           function(i_lambda) {
             logLik <- rep(NA_real_, nrow(data))
             for(k in seq_len(K))
               logLik[CV_folds == k] <- l_results_CV[[k]]$l_logLik[[i_lambda]]
             return(logLik)
           })
  ll_fits_CV <- lapply(l_results_CV, function(result_k) result_k$l_fits)
  
  return(list(l_fits_full = l_fits_full,
              logLik_CV = logLik_CV,
              ll_fits_CV = ll_fits_CV,
              l_filtering = l_filtering,
              CV_folds = CV_folds))
}


#' Simulation of realistic microbial observations
#'
#' @param n_sample Number of samples to simulate
#' @param n_feature Number of features to simulate
#' @param template Template for simulation. This can be either one of the prescribed datasets or fitting results from SparseDOSSA2_fit
#' @param new_features Should SparseDOSSA2 simulate new microbial features, or the same as the template dataset?
#' @param spike_metadata Should SparseDOSSA2 simulate associations with metadata?
#' @param spike_strength Effect size of association with metadata
#' @param control 
#'
#' @return
#' @export
#'
#' @examples
SparseDOSSA2_generate <- function(n_sample = NULL,
                                  n_feature = NULL,
                                  template = NULL,
                                  same_features = FALSE,
                                  spike_metadata = NULL,
                                  spike_strength = NULL,
                                  spike_feature = NULL,
                                  spike_discrete_bidirectional = TRUE,
                                  read_depth = 50000,
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
                                              n_feature = n_feature)
    rownames(featureParam_new) <- paste0("Feature", 1:n_feature)
  }
  
  # generate basis matrix
  if(verbose) message("Generating basis relative abundance matrix...")
  lMat_basis <- generate_basis(feature_param = featureParam_new,
                               n_sample = n_sample,
                               fit_C = fit_C)
  
  
  # if(!is.null(spike_metadata)) {
  #w
  # }
  
  if(verbose) message("Generating count matrix...")
  # generate read depth
  readDepth_new <- generate_readDepth(mu = fit_readCount["mu"],
                                      sigma2 = fit_readCount["sigma2"],
                                      n_sample = n_sample,
                                      read_depth = read_depth)
  # generate read counts
  mat_count <- generate_readCount(mat_p = lMat_basis$mat_basis,
                                  n_i = readDepth_new)
  
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
