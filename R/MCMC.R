EM <- function(data, controls = list(ncores = 6,
                                     lambda = 0.2,
                                     R = 1000,
                                     burnin = 0.1,
                                     threshold = 1e-3,
                                     maxit = 1000,
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
  l_easums <- list()
  ll_params <- list()
  while(TRUE) {
    i_iter <- i_iter + 1
    if(controls$verbose)
      print(i_iter)
    
    ## E step
    doParallel::registerDoParallel(cores = controls$ncores)
    # e_asums <- foreach::`%dopar%`(
    #   foreach::foreach(i_sample = seq_len(nrow(data)),
    #                    .combine='c'),
    #   {
    #     asum_samples <- mcmc_asum(x = data[i_sample, , drop = TRUE],
    #                               params = params,
    #                               R = controls$R)
    #     mean(vapply(
    #       asum_samples[-seq_len(round(controls$R * controls$burnin))],
    #       function(i_asum) i_asum[["asum"]],
    #       0.0))
    #   })
    e_asums <- foreach::`%dopar%`(
      foreach::foreach(i_sample = seq_len(nrow(data)),
                       .combine='c'),
      {
        asum_num <- integrate(vintegrand_num_asum,
                              lower = -20, upper = 20, ## FIXME
                              x = data[i_sample, , drop = TRUE],
                              params = params)
        asum_denom <- integrate(vintegrand_denom_asum,
                                lower = -20, upper = 20, ## FIXME
                                x = data[i_sample, , drop = TRUE],
                                params = params)
        asum_num$value / asum_denom$value
      })
    doParallel::stopImplicitCluster()
    l_easums[[i_iter]] <- e_asums
    
    ## M step
    a_data <- data * e_asums
    # a_data[is.na(e_asums), ] <- 
    #   data[is.na(e_asums), ] * 
    #   mean(e_asums, na.rm = TRUE) ## FIXME
    fit_sigmas <- get_sigmas(a_data, params$mu)
    fit_copulasso <- copulasso(data = a_data, 
                               lambda_list = controls$lambda,
                               K_CV = NULL, ## FIXME
                               ncores = controls$ncores
    )
    params_new <- list(pi0 = fit_marginals[, 1],
                       mu = fit_marginals[,2],
                       sigma = fit_sigmas,
                       Omega = fit_copulasso$fits[[1]])
    
    diff <- vapply(seq(from = 3, to = length(params_new)),
                   function(i_param) {
                     max(abs(params_new[[i_param]] - params[[i_param]]))
                   },
                   0.0)
    print(diff)
    
    if(all(diff < controls$threshold))
      break
    if(i_iter > controls$maxit)
      stop("Maximum iteration reached!")
    
    ll_params[[i_iter]] <- params_new
    params <- params_new
  }
}


mcmc_asum <- function(x,
                      params,
                      R = 100) {
  l_samples <- list()
  for(r in seq_len(R)) {
    if(r == 1)
      l_samples[[r]] <- 
        one_step_asum(asum = NULL, logLik = NULL, 
                      x = x, params = params)
    else 
      l_samples[[r]] <- 
        one_step_asum(asum = l_samples[[r - 1]]$asum, 
                      logLik = l_samples[[r - 1]]$logLik,
                      x = x, params = params)
  }
  
  return(l_samples)
}

one_step_asum <- function(asum = NULL, logLik = NULL, 
                          x, params) {
  
  if(is.null(asum)) {
    asum_star <- sum(x) ## FIXME?
    logLik <- -Inf
  }
  else 
    asum_star <- exp(log(asum) + rnorm(n = 1, sd = 1)) ## FIXME
  
  u_star <- a_to_u(a(x = x, asum = asum_star), 
                   pi0 = params$pi0, mu = params$mu, sigma = params$sigma)
  g_star <- u_to_g(u_star, a = a(x = x, asum = asum_star), 
                   mu = params$mu, sigma = params$sigma)
  logLik_star <- logLik_copula(g = g_star, asum = asum_star, x = x,
                               mu = params$mu, sigma = params$sigma, 
                               Omega = params$Omega)
  
  if(logLik_star == -Inf) ## FIXME
    prop <- 0
  else
    prop <- exp(logLik_star - logLik)
  
  if(runif(1) < prop)
    return(list(asum = asum_star,
                logLik = logLik_star,
                accept = TRUE))
  else
    return(list(asum = asum,
                logLik = logLik,
                accept = FALSE))
}


