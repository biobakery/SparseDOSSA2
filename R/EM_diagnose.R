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
  
  # Initialize using relative abundances
  fit_marginals <- get_marginals(data)
  
  l_results <- list()
  for(i_lambda in seq_along(lambdas)) {
    if(control$verbose) message("Fitting for lambda ", i_lambda)
    lambda <- lambdas[i_lambda]
    
    l_results[[i_lambda]] <- future::future({
      
      # Initialize using relative abundances
      fit_copulasso <- copulasso(data = data, 
                                 lambda_list = lambda,
                                 K_CV = NULL) ## FIXME
      params <- list(pi0 = fit_marginals[, 1],
                     mu = fit_marginals[,2],
                     sigma = fit_marginals[, 3],
                     Sigma = solve(fit_copulasso$fits[[1]]),
                     Omega = fit_copulasso$fits[[1]])
      
      i_iter <- 0
      converge <- FALSE
      ll_easums <- list()
      ll_params <- list()
      while(TRUE) {
        time <- Sys.time()
        if(i_iter + 1 > control$maxit) break
        i_iter <- i_iter + 1
        if(control$verbose)
          message("EM iteration ", i_iter)
        
        ## E step
        e_asums <- future.apply::future_vapply(
          seq_len(nrow(data)),
          function(i_sample) {
            i_time <- Sys.time()
            if(i_iter == 1) offset_a <- 1
            else offset_a <- ll_easums[[i_iter - 1]][i_sample, 1]
            num <- ea(x = data[i_sample, , drop = TRUE],
                      pi0 = params$pi0, mu = params$mu, 
                      sigma = params$sigma, Omega = params$Omega,
                      offset_a = offset_a,
                      control = c(control$control_numint, 
                                  list(only_value = FALSE,
                                       proper = FALSE)))
            eloga_num <- eloga(x = data[i_sample, , drop = TRUE],
                               pi0 = params$pi0, mu = params$mu, 
                               sigma = params$sigma, Omega = params$Omega,
                               offset_a = offset_a,
                               control = c(control$control_numint, 
                                           list(only_value = FALSE,
                                                proper = FALSE)))
            eloga2_num <- eloga2(x = data[i_sample, , drop = TRUE],
                                 pi0 = params$pi0, mu = params$mu, 
                                 sigma = params$sigma, Omega = params$Omega,
                                 offset_a = offset_a,
                                 control = c(control$control_numint, 
                                             list(only_value = FALSE,
                                                  proper = FALSE)))
            denom <- dx(x = data[i_sample, , drop = TRUE],
                        pi0 = params$pi0, mu = params$mu, 
                        sigma = params$sigma, Omega = params$Omega,
                        offset_a = offset_a,
                        control = c(control$control_numint, 
                                    list(only_value = FALSE,
                                         proper = FALSE)))
            l <- log_dx(x = data[i_sample, , drop = TRUE],
                        pi0 = params$pi0, mu = params$mu,
                        sigma = params$sigma, Omega = params$Omega,
                        offset_a = offset_a,
                        control = control$control_numint)
            return(c("mean" = num$integral / denom$integral,
                     "error" = abs(num$error / denom$integral) + 
                       abs(num$integral / (denom$integral)^2 * 
                             denom$error),
                     "l" = l,
                     "eloga" = eloga_num$integral / denom$integral,
                     "eloga2" = eloga2_num$integral / denom$integral,
                     "time" = Sys.time() - i_time))
          },
          rep(0.0, 6)
        ) %>% t()
        ll_easums[[i_iter]] <- e_asums
        
        ## M step
        a_data <- (data * e_asums[, 1])[!is.na(e_asums[, 1]), ] ## FIXME
        fit_sigmas <- get_sigmas(x = data, 
                                 eloga = e_asums[, "eloga"], 
                                 eloga2 = e_asums[, "eloga2"], 
                                 mu = fit_marginals[, 2])
        fit_copulasso <- copulasso(data = a_data, 
                                   lambda_list = lambda,
                                   K_CV = NULL) ## FIXME
        params_new <- list(pi0 = fit_marginals[, 1],
                           mu = fit_marginals[, 2],
                           sigma = fit_sigmas,
                           Sigma = solve(fit_copulasso$fits[[1]]),
                           Omega = fit_copulasso$fits[[1]])
        
        diff_abs <- vapply(c("sigma", "Sigma"), 
                           function(i_param)
                             get_diff(params_new[[i_param]], params[[i_param]], 
                                      denom_c = control$abs_tol, method = "abs"),
                           0.0)
        diff_rel <- vapply(c("sigma", "Sigma"), 
                           function(i_param)
                             get_diff(params_new[[i_param]], params[[i_param]], 
                                      denom_c = control$abs_tol, method = "rel"),
                           0.0)
        
        ll_params[[i_iter]] <- c(params_new,
                                 list(diff = c(diff_abs, diff_rel),
                                      l = mean(ll_easums[[i_iter]][, 3]),
                                      time = Sys.time() - time))
        params <- params_new
        
        if(!is.null(control$debug_dir)) {
          l_debug <- list(ll_easums = ll_easums, ll_params = ll_params)
          save(l_debug,
               file = paste0(control$debug_dir,"debug_lambda_", i_lambda, ".RData"))
        }
        
        if(max(diff_abs) < control$abs_tol & max(diff_rel) < control$rel_tol) {
          converge <- TRUE
          break
        }
      }
      
      list(lambda = lambda,
           fit = ll_params[[i_iter]],
           convergence = list(converge = converge,
                              n_iter = i_iter))
    })
  }
  
  l_results <- future::values(l_results)
  return(l_results)
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
  
  CV_folds <- make_CVfolds(n = nrow(data), K = K)
  
  l_results_CV <- list()
  for(k in seq_len(K)) {
    if(control$verbose) message("Performing CV K=", k)
    l_results_CV[[k]] <- future::future({
      data_training <- data[CV_folds != k, ]
      data_testing <- data[CV_folds == k, ]
      control_tmp <- control
      control_tmp$debug_dir <- paste0(control$debug_dir, "K", k, "/")
      
      l_fits <- EM_diagnose(data = data_training,
                            lambdas = lambdas,
                            control = control_tmp)
      
      l_ll <- future.apply::future_lapply(
        seq_along(lambdas),
        function(i_lambda) {
          params <- l_fits[[i_lambda]]$fit
          future.apply::future_sapply(
            seq_len(nrow(data_testing)),
            function(i_sample) {
              ll <- log_dx(x = data_testing[i_sample, , drop = TRUE],
                           pi0 = params$pi0, mu = params$mu,
                           sigma = params$sigma, Omega = params$Omega,
                           offset_a = 1,
                           control = control$control_numint)
            })
        })
      
      for(i_lambda in seq_along(l_fits)) {
        l_fits[[i_lambda]]$ll_CV <- l_ll[[i_lambda]]
      }
      l_fits
    })
  }
  
  if(control$verbose) message("Aggregating likelihood")
  ll_CV <- sapply(seq_len(lambdas),
                  function(i_lambda) {
                    ll <- rep(NA_real_, nrow(data))
                    for(k in seq_len(K))
                      ll[CV_folds == k] <- l_results_CV[[k]][[i_lambda]]$ll_CV
                    return(ll)
                  })
  
  if(control$verbose) message("Performing overall fit")
  control_tmp <- control
  control_tmp$debug_dir <- paste0(control_tmp$debug_dir, "K0/")
  l_results <- future::future({
    EM_diagnose(data = data_training,
                lambdas = lambdas,
                control = control_tmp)
  })
  
  return(list(full = l_results,
              ll_CV = ll_CV,
              CV = l_results_CV))
}

control_EM <- function(maxit = 1000,
                       rel_tol = 1e-3,
                       abs_tol = 1e-3,
                       control_numint = list(subdivisions = 10000,
                                             limit_max = 50,
                                             limit_min = 1e-10,
                                             step_size = 2,
                                             rel_tol = 1e-6,
                                             abs_tol = 0),
                       verbose = FALSE,
                       debug_dir = NULL) {
  list(maxit = maxit,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       control_numint = control_numint,
       verbose = verbose,
       debug_dir = debug_dir)
}