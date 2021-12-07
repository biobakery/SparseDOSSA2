EM <- function(data, 
               lambda,
               control = list()) {
  # initialization, filtering
  if(is.null(lambda) & any(lambda <= 0))
    stop("lambda must be a positive value!")
  control <- do.call(control_fit, control)
  debug_copulasso_file <- NULL
  if(!is.null(control$debug_dir)) {
    dir.create(control$debug_dir, recursive = TRUE)
    control$debug_dir <- normalizePath(control$debug_dir)
  }
    
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  # initialize EM using relative abundances
  feature_param <- fit_featureParam(data)
  fit_copulasso <- copulasso(
    data = data, 
    marginals = feature_param,
    lambda = lambda,
    control = list(debug_dir = control$debug_dir)
  )
  feature_param[, "mu"] <- feature_param[, "mu"] - mean(feature_param[, "mu"])
  params <- list(pi0 = feature_param[ ,"pi0"],
                 mu = feature_param[, "mu"],
                 sigma = feature_param[, "sigma"],
                 Sigma = fit_copulasso$Sigma,
                 Omega = fit_copulasso$Omega,
                 Corr_star = fit_copulasso$Corr_star,
                 diff = rep(NA_real_, 4),
                 logLik = -Inf,
                 time = Sys.time())
  if(fit_copulasso$copulasso_code != 0) {
    warning("Missing values in Omega estimation! (lambda to small?)")
    return(list(lambda = lambda,
                fit = params,
                convergence = list(converge = FALSE,
                                   converge_code = 4,
                                   n_iter = i_iter)))
  }
  
  # EM algorithim    
  i_iter <- 0
  converge <- FALSE
  ll_easums <- list()
  ll_params <- list()
  while(TRUE) {
    i_iter <- i_iter + 1
    if(control$verbose)
      message("EM iteration ", i_iter)
    
    # E step
    e_asums <- future.apply::future_vapply(
      seq_len(nrow(data)),
      function(i_sample)
        get_es(x = data[i_sample, , drop = TRUE],
               pi0 = params$pi0, mu = params$mu, sigma = params$sigma, 
               Omega = params$Omega, Sigma = params$Sigma,
               control = control$control_numint),
      rep(0.0, 10)
    ) %>% t()
    if(any(is.na(e_asums[, c("ea", "dx", "eloga", "eloga2")]))) {
      warning("Numeric integration in E step returned NAs!")
      converge_code <- 3
      break
    }
    if(any(e_asums[, "eloga"]^2 >= e_asums[, "eloga2"])) {
      warning("Numeric integration in E step gave bad expectation values!")
      converge_code <- 2
      break
    }
    ll_easums[[i_iter]] <- e_asums
    
    ## M step
    fit_copulasso <- copulasso(data = data * exp(e_asums[, "eloga"]), 
                               marginals = feature_param,
                               lambda = lambda)
    if(fit_copulasso$copulasso_code != 0) {
      warning("Missing values in Omega estimation! (lambda to small?)")
      converge_code <- 4
      break
    }
    sigmahat <- get_sigmahat(data = data,
                             eloga = e_asums[, "eloga"],
                             eloga2 = e_asums[, "eloga2"])
    feature_param[, "sigma"] <- sigmahat
    params_new <- list(pi0 = feature_param[ ,"pi0"],
                       mu = feature_param[ ,"mu"],
                       sigma = feature_param[ ,"sigma"],
                       Sigma = fit_copulasso$Sigma,
                       Omega = fit_copulasso$Omega,
                       Corr_star = fit_copulasso$Corr_star,
                       logLik = mean(e_asums[, "logLik"]))
    diff_abs <- get_diff(params_new[["logLik"]], params[["logLik"]], 
                         denom_c = control$abs_tol, method = "abs")
    diff_rel <- get_diff(params_new[["logLik"]], params[["logLik"]], 
                         denom_c = control$abs_tol, method = "rel")
    params <- c(params_new,
                list(diff = c(diff_abs, diff_rel),
                     time = Sys.time()))
    ll_params[[i_iter]] <- params
    
    if(!is.null(control$debug_dir)) {
      l_debug <- list(ll_easums = ll_easums, 
                      ll_params = ll_params, 
                      l_filtering = l_filtering)
      save(l_debug,
           file = paste0(control$debug_dir,
                         "/debug_EM.RData"))
    }
    
    if(max(diff_abs) < control$abs_tol & max(diff_rel) < control$rel_tol) {
      converge <- TRUE
      converge_code <- 0
      break
    }
    if(i_iter + 1 > control$maxit) {
      warning("Maximum EM iteration reached!")
      converge_code <- 1
      break
    }
  }
  
  return(list(fit = params,
              lambda = lambda,
              convergence = list(converge = converge,
                                 converge_code = converge_code,
                                 n_iter = i_iter),
              l_filtering = l_filtering))
}

EM_CV <- function(data,
                  lambdas,
                  K = 5,
                  control = list()) {
  # initialization, filtering
  control <- do.call(control_fit, control)
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  if(is.null(lambdas))
    lambdas <- 10^seq(-2, 0, length.out = 5)
  if(any(lambdas <= 0))
    stop("lambdas must be positive values!")
  if(!is.null(control$debug_dir)) {
    dir.create(control$debug_dir, recursive = TRUE)
    control$debug_dir <- normalizePath(control$debug_dir)
  }

  if(control$verbose) message("Performing full data fit...")
  l_fits_full <- future::future({
    future.apply::future_lapply(
      seq_along(lambdas),
      function(i_lambda) {
        control_tmp <- control
        control_tmp$verbose <- FALSE
        if(!is.null(control$debug_dir))
          control_tmp$debug_dir <- paste0(control$debug_dir,
                                          "/K_0/lambda_",
                                          i_lambda)
        EM(data = data,
           lambda = lambdas[i_lambda],
           control = control_tmp)
      })
  })
  l_fits_full <- future::value(l_fits_full)

  # CV fits
  CV_folds <- make_CVfolds(n = nrow(data), K = K)
  if(!is.null(control$debug_dir))
    save(CV_folds, file = paste0(control$debug_dir, "/CV_folds.RData"))
  ll_results_CV <- list()
  for(k in seq_len(K)) {
    if(control$verbose) 
      message("Performing CV k=", k)
    ll_results_CV[[k]] <- future::future({
      data_training <- data[CV_folds != k, , drop = FALSE]
      data_testing <- data[CV_folds == k, , drop = FALSE]
      
      l_fits_k <- 
        future.apply::future_lapply(
          seq_along(lambdas),
          function(i_lambda) {
            control_tmp <- control
            control_tmp$verbose <- FALSE
            if(!is.null(control$debug_dir))
              control_tmp$debug_dir <- paste0(control$debug_dir,
                                              "/K_", k,
                                              "/lambda_", i_lambda)
            EM(data = data_training,
               lambda = lambdas[i_lambda],
               control = control_tmp)
          })
      
      # Fill in parameters estimates for features not present in training data
      for(i_lambda in seq_along(lambdas))
        l_fits_k[[i_lambda]]$fit <-
        fill_estimates_CV(l_fits_k[[i_lambda]]$fit,
                          l_fits_full[[i_lambda]]$fit,
                          l_fits_k[[i_lambda]]$l_filtering$ind_feature)

      # Calculate ll in testing data
      l_logLik <- future.apply::future_lapply(
        seq_along(lambdas),
        function(i_lambda) {
          params <- l_fits_k[[i_lambda]]$fit
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

      list(l_fits = l_fits_k,
           l_logLik = l_logLik)
    })
  }
  ll_results_CV <- future::value(ll_results_CV)

  # Aggregate CV results
  logLik_CV <-
    sapply(seq_along(lambdas),
           function(i_lambda) {
             logLik <- rep(NA_real_, nrow(data))
             for(k in seq_len(K))
               logLik[CV_folds == k] <- ll_results_CV[[k]]$l_logLik[[i_lambda]]
             return(logLik)
           })
  ll_fits_CV <- lapply(ll_results_CV, function(results_k) results_k$l_fits)

  return(list(fit = 
                l_fits_full[[
                  order(
                    apply(-logLik_CV, 
                          2, 
                          function(x)
                            mean(setdiff(x, Inf))))[1]]]$fit,
              lambdas = lambdas,
              logLik_CV = logLik_CV,
              l_fits_full = l_fits_full,
              CV_folds = CV_folds,
              ll_fits_CV = ll_fits_CV,
              l_filtering = l_filtering
  ))
}
