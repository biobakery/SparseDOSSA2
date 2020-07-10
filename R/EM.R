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
EM <- function(data, 
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
          rep(0.0, 10)
        ) %>% t()
        if(any(is.na(e_asums[, c("ea", "dx", "eloga", "eloga2")]))) {
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
