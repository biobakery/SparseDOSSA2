#' Title
#'
#' @param data 
#' @param params_init 
#' @param control 
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
EM_diagnose <- function(data, 
                        lambdas,
                        control = list()) {
  control <- do.call(control_EM, control)
  if(!is.null(control$debug_dir))
    dir.create(control$debug_dir)
  
  # Filtering out features/samples that are all zero
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  
  # Initialize using relative abundances
  fit_marginals <- get_marginals(data)
  
  l_fits <- list()
  for(i_lambda in seq_along(lambdas)) {
    if(control$verbose) message("Fitting for lambda ", i_lambda)
    lambda <- lambdas[i_lambda]
    
    l_fits[[i_lambda]] <- future::future({
      
      # Initialize using relative abundances
      fit_copulasso <- copulasso(
        data = data, 
        marginals = fit_marginals,
        lambda = lambda,
        control = 
          c(control$control_copulasso,
            list(debug_file = 
                   paste0(control$debug_dir,
                          "debug_glasso_lambda_", 
                          i_lambda, ".RData")))) ## FIXME
      params <- list(pi0 = fit_marginals$pi0,
                     mu = fit_marginals$mu,
                     sigma = fit_marginals$sigma,
                     Sigma = fit_copulasso$Sigma,
                     Omega = fit_copulasso$Omega,
                     Corr_star = fit_copulasso$Corr_star,
                     diff = rep(NA_real_, 4),
                     logLik = NA_real_,
                     time = Sys.time())
      
      i_iter <- 0
      converge <- FALSE
      if(fit_copulasso$copulasso_code != 0) {
        warning("Missing values in Omega estimation! (lambda to small?)")
        return(list(lambda = lambda,
                    fit = params,
                    convergence = list(converge = converge,
                                       converge_code = 4,
                                       n_iter = i_iter)))
      }
        
      ll_easums <- list()
      ll_params <- list()
      while(TRUE) {
        i_iter <- i_iter + 1
        if(control$verbose)
          message("EM iteration ", i_iter)
        
        ## E step
        e_asums <- future.apply::future_vapply(
          seq_len(nrow(data)),
          function(i_sample)
            get_es(x = data[i_sample, , drop = TRUE],
                   pi0 = params$pi0, mu = params$mu, sigma = params$sigma, 
                   Omega = params$Omega, Sigma = params$Sigma,
                   control = control$control_numint),
          rep(0.0, 9)
        ) %>% t()
        if(any(is.na(e_asums[, c("ea", "logLik", "eloga", "eloga2")]))) {
          warning("Numeric integration in E step returned NAs!")
          converge_code <- 3
          break
        }
        if(any(e_asums[, "eloga"]^2 > e_asums[, "eloga2"])) {
          warning("Numeric integration in E step gave bad expectation values!")
          converge_code <- 2
          break
        }
        ll_easums[[i_iter]] <- e_asums
        
        ## M step
        a_data <- data * e_asums[, "ea"] ## FIXME
        fit_sigmas <- get_sigmas(x = data, 
                                 eloga = e_asums[, "eloga"], 
                                 eloga2 = e_asums[, "eloga2"], 
                                 mu = fit_marginals$mu)
        fit_marginals$sigma <- fit_sigmas
        fit_copulasso <- copulasso(data = a_data, 
                                   marginals = fit_marginals,
                                   lambda = lambda,
                                   control = control$control_copulasso)
        if(fit_copulasso$copulasso_code != 0) {
          warning("Missing values in Omega estimation! (lambda to small?)")
          converge_code <- 4
          break
        }
        params_new <- list(pi0 = fit_marginals$pi0,
                           mu = fit_marginals$mu,
                           sigma = fit_sigmas,
                           Sigma = fit_copulasso$Sigma,
                           Omega = fit_copulasso$Omega,
                           Corr_star = fit_copulasso$Corr_star)
        diff_abs <- vapply(
          c("sigma", "Corr_star"), 
          function(i_param)
            get_diff(params_new[[i_param]], params[[i_param]], 
                     denom_c = control$abs_tol, method = "abs"),
          0.0)
        diff_rel <- vapply(
          c("sigma", "Corr_star"), 
          function(i_param)
            get_diff(params_new[[i_param]], params[[i_param]], 
                     denom_c = control$abs_tol, method = "rel"),
          0.0)
        params <- c(params_new,
                    list(diff = c(diff_abs, diff_rel),
                         logLik = mean(e_asums[, "logLik"]),
                         time = Sys.time()))
        ll_params[[i_iter]] <- params
        
        if(!is.null(control$debug_dir)) {
          l_debug <- list(ll_easums = ll_easums, 
                          ll_params = ll_params, 
                          l_filtering = l_filtering)
          save(l_debug,
               file = paste0(control$debug_dir,
                             "debug_lambda_", i_lambda, ".RData"))
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
      
      return(list(lambda = lambda,
                  fit = params,
                  convergence = list(converge = converge,
                                     converge_code = converge_code,
                                     n_iter = i_iter)))
    })
  }
  
  l_fits <- future::values(l_fits)
  return(list(l_fits = l_fits,
              l_filtering = l_filtering))
}


#' Title
#'
#' @param data 
#' @param params_init 
#' @param control 
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
EM_diagnose_CV <- function(data, 
                           lambdas,
                           K,
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
    EM_diagnose(data = data,
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
      
      result_k <- EM_diagnose(data = data_training,
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
              logLik <- dx(x = data_testing[i_sample, , drop = TRUE],
                           pi0 = params$pi0, mu = params$mu,
                           sigma = params$sigma, Omega = params$Omega,
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
  logLik_CV <- sapply(seq_along(lambdas),
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

control_EM <- function(maxit = 100,
                       rel_tol = 1,
                       abs_tol = 1e-2,
                       control_numint = list(),
                       control_copulasso = list(),
                       verbose = FALSE,
                       debug_dir = NULL) {
  list(maxit = maxit,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       control_numint = control_numint,
       control_copulasso = control_copulasso,
       verbose = verbose,
       debug_dir = debug_dir)
}