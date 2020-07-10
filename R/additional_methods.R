mvn <- function(data, 
                lambdas) {
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  data_pseudo <- t(apply(data + min(setdiff(data, 0)) / 2,
                         1,
                         function(x) x / sum(x)))
  data_transformed <- log(data_pseudo) - log(1 - data_pseudo)
  
  l_fits <- list()
  for(i_k in seq_along(lambdas)) {
    lambda <- lambdas[i_k]
    
    l_fits[[i_k]] <- future::future({
      
      # marginal parameters
      mu <- apply(data_transformed, 2, mean)
      sigma <- apply(data_transformed, 2, sd)
      
      # glasso fit
      S <- cor(data_transformed, method = "pearson")
      S_conditioned <- condition_ridge(S,
                                       lambda = 1e-6,
                                       method = "ridge1")
      Omega <- huge::huge.glasso(x = S_conditioned,
                                 lambda = lambda,
                                 verbose = FALSE)$icov[[1]]
      
      Omega <- enforce_symm(Omega, method = "svd")
      Omega <- enforce_corr(Omega)
      
      Sigma <- threshold_matrix(solve(Omega), 
                                threshold_zero = 1e-16)
      
      
      return(list(lambda = lambda,
                  fit = list(mu = mu, 
                             sigma = sigma,
                             Sigma = Sigma)))
    })
  }
  
  l_fits <- future::values(l_fits)
  return(list(l_fits = l_fits,
              l_filtering = l_filtering))
}

mvn_CV <- function(data, 
                   lambdas,
                   K) {
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  l_fits_full <- future::future({
    mvn(data = data,
        lambdas = lambdas)$l_fits
  })
  l_fits_full <- future::value(l_fits_full)
  
  # CV fits
  CV_folds <- make_CVfolds(n = nrow(data), K = K)
  l_results_CV <- list()
  for(k in seq_len(K)) {
    l_results_CV[[k]] <- future::future({
      data_training <- data[CV_folds != k, ]
      data_testing <- data[CV_folds == k, ]
      data_testing_pseudo <- 
        t(apply(data_testing + min(setdiff(data_testing, 0)) / 2,
                1,
                function(x) x / sum(x)))
      data_testing_transformed <- 
        log(data_testing_pseudo) - 
        log(1 - data_testing_pseudo)
      
      result_k <- mvn(data = data_training,
                      lambdas = lambdas)
      
      # Fill in parameters estimates for features not present in training data
      for(i_k in seq_along(lambdas))
        result_k$l_fits[[i_k]]$fit <- 
        fill_estimates_CV_mvn(result_k$l_fits[[i_k]]$fit,
                              l_fits_full[[i_k]]$fit,
                              result_k$l_filtering$ind_feature)
      
      # Calculate ll in testing data
      l_logLik <- future.apply::future_lapply(
        seq_along(lambdas),
        function(i_k) {
          params <- result_k$l_fits[[i_k]]$fit
          future.apply::future_sapply(
            seq_len(nrow(data_testing)),
            function(i_sample) {
              logLik <- 
                mvtnorm::dmvnorm(
                  x = data_testing_pseudo[i_sample, , drop = TRUE],
                  mean = params$mu, 
                  sigma = diag(sqrt(params$sigma)) %*%
                    params$Sigma %*%
                    diag(sqrt(params$sigma)),
                  log = TRUE)
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
           function(i_k) {
             logLik <- rep(NA_real_, nrow(data))
             for(k in seq_len(K))
               logLik[CV_folds == k] <- l_results_CV[[k]]$l_logLik[[i_k]]
             return(logLik)
           })
  ll_fits_CV <- lapply(l_results_CV, function(result_k) result_k$l_fits)
  
  return(list(l_fits_full = l_fits_full,
              logLik_CV = logLik_CV,
              ll_fits_CV = ll_fits_CV,
              l_filtering = l_filtering,
              CV_folds = CV_folds))
}

fill_estimates_CV_mvn <- function(params_CV, params_full, ind_feature) {
  params_return <- params_CV
  
  params_return[c("mu", "sigma")] <- params_full[c("mu", "sigma")]
  for(param in c("mu", "sigma"))
    params_return[[param]][ind_feature] <- params_CV[[param]]
  
  params_return$Sigma <- 
    diag(rep(1, length(ind_feature)))
  params_return[["Sigma"]][ind_feature, ind_feature] <- params_CV[["Sigma"]]
  
  return(params_return)
}

dmm <- function(data, 
                ks) {
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  l_fits <- list()
  for(i_k in seq_along(ks)) {
    k <- ks[i_k]
    
    l_fits[[i_k]] <- future::future({
      
      # marginal parameters
      dmm_fit <- DirichletMultinomial::dmn(count = data,
                                           k = k, 
                                           verbose = FALSE)
      
      return(list(k = k,
                  fit = list(wt = DirichletMultinomial::mixturewt(dmm_fit)[, 1], 
                             alpha = t(DirichletMultinomial::fitted(dmm_fit)))))
    })
  }
  
  l_fits <- future::values(l_fits)
  return(list(l_fits = l_fits,
              l_filtering = l_filtering))
}

dmm_CV <- function(data, 
                   ks,
                   K) {
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  l_fits_full <- future::future({
    dmm(data = data,
        ks = ks)$l_fits
  })
  l_fits_full <- future::value(l_fits_full)
  
  # CV fits
  CV_folds <- make_CVfolds(n = nrow(data), K = K)
  l_results_CV <- list()
  for(k in seq_len(K)) {
    l_results_CV[[k]] <- future::future({
      data_training <- data[CV_folds != k, ]
      data_testing <- data[CV_folds == k, ]
      
      result_k <- dmm(data = data_training,
                      ks = ks)
      
      # Fill in parameters estimates for features not present in training data
      for(i_k in seq_along(ks))
        result_k$l_fits[[i_k]]$fit <- 
        fill_estimates_CV_dmm(result_k$l_fits[[i_k]]$fit,
                              l_fits_full[[i_k]]$fit,
                              result_k$l_filtering$ind_feature)
      
      # Calculate ll in testing data
      l_logLik <- future.apply::future_lapply(
        seq_along(ks),
        function(i_k) {
          params <- result_k$l_fits[[i_k]]$fit
          future.apply::future_sapply(
            seq_len(nrow(data_testing)),
            function(i_sample) {
              logLik <- 
                ddmm(
                  x = data_testing[i_sample, , drop = TRUE],
                  wt = params$wt, 
                  alpha = params$alpha,
                  log = TRUE)
            })
        })
      
      list(l_fits = result_k$l_fits,
           l_logLik = l_logLik)
    })
  }
  l_results_CV <- future::values(l_results_CV)
  
  # Aggregate CV results
  logLik_CV <- 
    sapply(seq_along(ks),
           function(i_k) {
             logLik <- rep(NA_real_, nrow(data))
             for(k in seq_len(K))
               logLik[CV_folds == k] <- l_results_CV[[k]]$l_logLik[[i_k]]
             return(logLik)
           })
  ll_fits_CV <- lapply(l_results_CV, function(result_k) result_k$l_fits)
  
  return(list(l_fits_full = l_fits_full,
              logLik_CV = logLik_CV,
              ll_fits_CV = ll_fits_CV,
              l_filtering = l_filtering,
              CV_folds = CV_folds))
}

fill_estimates_CV_dmm <- function(params_CV, params_full, ind_feature) {
  params_return <- params_CV
  
  params_return[c("wt", "alpha")] <- params_full[c("wt", "alpha")]
  params_return[["alpha"]][, ind_feature] <- params_CV[["alpha"]]
  
  return(params_return)
}

ddmm <- function(x, wt, alpha, log = FALSE) {
  lik <- vapply(seq_along(wt),
                function(i_comp)
                  wt[i_comp] *
                  extraDistr::ddirmnom(x = x,
                                       size = sum(x),
                                       alpha = alpha[i_comp, ]),
                0.0)
  lik <- sum(lik)
  
  if(log)
    return(log(lik))
  return(lik)
}

ddm <- function(x, alpha) {
  factorial(sum(x)) * 
    gamma(sum(alpha)) /
    gamma(sum(x) + sum(alpha)) *
    prod(gamma(x + alpha) /
           factorial(x) /
           gamma(alpha))
}
