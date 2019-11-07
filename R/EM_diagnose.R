EM_diagnose <- function(data, controls = list(ncores = 6,
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
                             lambda_list = controls$lambda,
                             K_CV = NULL, ## FIXME
                             ncores = controls$ncores)
  params <- list(pi0 = fit_marginals[, 1],
                 mu = fit_marginals[,2],
                 sigma = fit_marginals[, 3],
                 Omega = fit_copulasso$fits[[1]])
  
  i_iter <- 0
  ll_easums <- list()
  ll_params <- list()
  while(TRUE) {
    i_iter <- i_iter + 1
    if(i_iter > controls$maxit) break
    if(controls$verbose)
      print(i_iter)
    
    ## E step
    doParallel::registerDoParallel(cores = controls$ncores)
    if(controls$method == "mcmc") {
      e_asums <- foreach::`%dopar%`(
        foreach::foreach(i_sample = seq_len(nrow(data)),
                         .combine='rbind'),
        {
          asum_samples <- mcmc_asum(x = data[i_sample, , drop = TRUE],
                                    params = params,
                                    R = controls$R)
          asum_samples <- vapply(
            asum_samples[-seq_len(round(controls$R * controls$burnin))],
            function(i_asum) i_asum[["asum"]],
            0.0)
          return(c("mean" = mean(asum_samples), 
                   "error" = sd(log(asum_samples)) / 
                     sqrt(controls$R*(1 - controls$burnin))))
        })
    }
    if(controls$method == "numint") {
      e_asums <- foreach::`%dopar%`(
        foreach::foreach(i_sample = seq_len(nrow(data)),
                         .combine='rbind'),
        {
          asum_num <- integrate(vintegrand_num_asum,
                                subdivisions = controls$subdivisions,
                                lower = -20, upper = 20, ## FIXME
                                x = data[i_sample, , drop = TRUE],
                                params = params)
          asum_denom <- integrate(vintegrand_denom_asum,
                                  subdivisions = controls$subdivisions,
                                  lower = -20, upper = 20, ## FIXME
                                  x = data[i_sample, , drop = TRUE],
                                  params = params)
          return(c("mean" = asum_num$value / asum_denom$value,
                   "error" = abs(asum_num$abs.error / asum_denom$value) + 
                     abs(asum_num$value / (asum_denom$value)^2 * 
                           asum_denom$abs.error)))
        })
    }
    doParallel::stopImplicitCluster()
    ll_easums[[i_iter]] <- e_asums
    
    ## M step
    a_data <- data * e_asums[, 1]
    fit_sigmas <- get_sigmas(a_data, params$mu)
    fit_copulasso <- copulasso(data = a_data, 
                               lambda_list = controls$lambda,
                               K_CV = NULL, ## FIXME
                               ncores = controls$ncores)
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
