EM_diagnose <- function(data, 
                        control = list()) {
  control <- do.call("control_EM", control)
  
  # Initialize using relative abundances
  fit_marginals <- get_marginals(data)
  fit_copulasso <- copulasso(data = data, 
                             lambda_list = control$lambda,
                             K_CV = NULL) ## FIXME
  params <- list(pi0 = fit_marginals[, 1],
                 mu = fit_marginals[,2],
                 sigma = fit_marginals[, 3],
                 Omega = fit_copulasso$fits[[1]])
  
  i_iter <- 0
  converge <- FALSE
  ll_easums <- list()
  ll_params <- list()
  while(TRUE) {
    i_iter <- i_iter + 1
    if(i_iter > control$maxit) break
    if(control$verbose)
      print(i_iter)
    
    ## E step
    if(control$method == "mcmc") {
      e_asums <- foreach::`%dopar%`(
        foreach::foreach(i_sample = seq_len(nrow(data)),
                         .combine='rbind'),
        {
          asum_samples <- mcmc_asum(x = data[i_sample, , drop = TRUE],
                                    params = params,
                                    R = control$control_mcmc$R)
          asum_samples <- vapply(
            asum_samples[-seq_len(round(control$control_mcmc$R * control$control_mcmc$burnin))],
            function(i_asum) i_asum[["asum"]],
            0.0)
          return(c("mean" = mean(asum_samples), 
                   "error" = sd(log(asum_samples)) / 
                     sqrt(control$control_mcmc$R*(1 - control$control_mcmc$burnin))))
        })
    }
    if(control$method == "numint") {
      e_asums <- future.apply::future_vapply(
        seq_len(nrow(data)),
        function(i_sample) {
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
          return(c("mean" = num$value / denom$value,
                   "error" = abs(num$abs.error / denom$value) + 
                     abs(num$value / (denom$value)^2 * 
                           denom$abs.error),
                   "l" = l,
                   "eloga" = eloga_num$value / denom$value,
                   "eloga2" = eloga2_num$value / denom$value))
        },
        rep(0.0, 5)
      ) %>% t()
    }
    ll_easums[[i_iter]] <- e_asums
    
    ## M step
    a_data <- (data * e_asums[, 1])[!is.na(e_asums[, 1]), ] ## FIXME
    fit_sigmas <- get_sigmas(x = data, 
                             eloga = e_asums[, "eloga"], 
                             eloga2 = e_asums[, "eloga2"], 
                             mu = fit_marginals[, 2])
    fit_copulasso <- copulasso(data = a_data, 
                               lambda_list = control$lambda,
                               K_CV = NULL) ## FIXME
    params_new <- list(pi0 = fit_marginals[, 1],
                       mu = fit_marginals[, 2],
                       sigma = fit_sigmas,
                       Omega = fit_copulasso$fits[[1]])
    
    diff_abs <- vapply(c(3, 4), 
                       function(i_param)
                         get_diff(params_new, params, 
                                  denom_c = control$abs_tol, method = "abs"),
                       0.0)
    diff_rel <- vapply(c(3, 4), 
                       function(i_param)
                         get_diff(params_new, params, 
                                  denom_c = control$abs_tol, method = "rel"),
                       0.0)
    
    ll_params[[i_iter]] <- c(params_new,
                             list(diff = c(diff_abs, diff_rel),
                                  l = sum(ll_easums[[i_iter]][, 3])))
    params <- params_new
    
    if(max(diff_abs) < control$abs_tol & max(diff_rel) < control$rel_tol) {
      converge <- TRUE
      break
    }
  }
  
  return(list(ll_easums = ll_easums, ll_params = ll_params,
              converge = list(converge = converge,
                              n_iter = i_iter)))
}

control_EM <- function(lambda = 0.2,
                       maxit = 30,
                       method = "numint",
                       rel_tol = 5e-2,
                       abs_tol = 1e-4,
                       control_mcmc = list(R = 10000,
                                           burnin = 0.1),
                       control_numint = list(subdivisions = 10000,
                                             limit_max = 50,
                                             limit_min = 1e-10,
                                             step_size = 2,
                                             rel.tol = 1e-6,
                                             abs.tol = 1e-6),
                       verbose = FALSE) {
  list(lambda = lambda,
       maxit = maxit,
       method = method,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       control_mcmc = control_mcmc,
       control_numint = control_numint,
       verbose = verbose)
}