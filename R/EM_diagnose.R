EM_diagnose <- function(data, 
                        control = list(ncores = 6,
                                       lambda = 0.2,
                                       R = 10000,
                                       subdivisions = 10000,
                                       burnin = 0.1,
                                       maxit = 30,
                                       method = "mcmc",
                                       verbose = TRUE)) {
  
  # Initialize using relative abundances
  fit_marginals <- get_marginals(data)
  fit_copulasso <- copulasso(data = data, 
                             lambda_list = control$lambda,
                             K_CV = NULL, ## FIXME
                             ncores = control$ncores)
  params <- list(pi0 = fit_marginals[, 1],
                 mu = fit_marginals[,2],
                 sigma = fit_marginals[, 3],
                 Omega = fit_copulasso$fits[[1]])
  
  i_iter <- 0
  ll_easums <- list()
  ll_params <- list()
  while(TRUE) {
    i_iter <- i_iter + 1
    if(i_iter > control$maxit) break
    if(control$verbose)
      print(i_iter)
    
    ## E step
    doParallel::registerDoParallel(cores = control$ncores)
    if(control$method == "mcmc") {
      e_asums <- foreach::`%dopar%`(
        foreach::foreach(i_sample = seq_len(nrow(data)),
                         .combine='rbind'),
        {
          asum_samples <- mcmc_asum(x = data[i_sample, , drop = TRUE],
                                    params = params,
                                    R = control$R)
          asum_samples <- vapply(
            asum_samples[-seq_len(round(control$R * control$burnin))],
            function(i_asum) i_asum[["asum"]],
            0.0)
          return(c("mean" = mean(asum_samples), 
                   "error" = sd(log(asum_samples)) / 
                     sqrt(control$R*(1 - control$burnin))))
        })
    }
    if(control$method == "numint") {
      e_asums <- foreach::`%dopar%`(
        foreach::foreach(i_sample = seq_len(nrow(data)),
                         .combine='rbind'),
        {
          num <- ea(x = data[i_sample, , drop = TRUE],
                    pi0 = params$pi0, mu = params$mu, 
                    sigma = params$sigma, Omega = params$Omega,
                    control = list(subdivisions = control$subdivisions,
                                   only_value = FALSE,
                                   proper = FALSE))
          denom <- dx(x = data[i_sample, , drop = TRUE],
                      pi0 = params$pi0, mu = params$mu, 
                      sigma = params$sigma, Omega = params$Omega,
                      control = list(subdivisions = control$subdivisions,
                                     only_value = FALSE,
                                     proper = FALSE))
          return(c("mean" = num$value / denom$value,
                   "error" = abs(num$abs.error / denom$value) + 
                     abs(num$value / (denom$value)^2 * 
                           denom$abs.error)))
        })
    }
    doParallel::stopImplicitCluster()
    ll_easums[[i_iter]] <- e_asums
    
    ## M step
    a_data <- (data * e_asums[, 1])[!is.na(e_asums[, 1]), ] ## FIXME
    fit_sigmas <- get_sigmas(a_data, params$mu)
    fit_copulasso <- copulasso(data = a_data, 
                               lambda_list = control$lambda,
                               K_CV = NULL, ## FIXME
                               ncores = control$ncores)
    params_new <- list(pi0 = fit_marginals[, 1],
                       mu = fit_marginals[, 2],
                       sigma = fit_sigmas,
                       Omega = fit_copulasso$fits[[1]])
    
    diff <- vapply(seq(from = 3, to = length(params_new)),
                   function(i_param) {
                     max(abs(params_new[[i_param]] - params[[i_param]]))
                   },
                   0.0)
    
    ll_params[[i_iter]] <- c(params_new, list("diff" = diff))
    params <- params_new
  }
  
  return(list(ll_easums = ll_easums, ll_params = ll_params))
}
