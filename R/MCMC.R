mcmc_a <- function(x,
                   params,
                   R = 100) {
  p <- length(x)
  
  a <- 1
  u <- a_to_u(a * x, pi0 = params$pi0, mu = params$mu, sigma = params$sigma)
  g <- qnorm(u)
  
  l_samples <- list()
  l_samples[[1]] <- 
    list(a = a,
         loglik = 
           log_dmvnorm(g = g, Omega = params$Omega) +
           log(a) * p,
         accept = TRUE)
  
  for(r in seq_len(R)) {
    l_samples[[r + 1]] <- 
      one_step_a(
        x = x,
        a = l_samples[[r]]$a,
        params = params,
        loglik = l_samples[[r]]$loglik)
  }
  
  return(l_samples)
}

one_step_a <- function(x, a, loglik, params) {
  p <- length(x)
  
  a_star <- exp(log(a) + rnorm(n = 1, sd = 0.1))
  u_star <- a_to_u(a_star * x, 
                   pi0 = params$pi0, mu = params$mu, sigma = params$sigma)
  g_star <- qnorm(u_star)
  loglik_star <- log_dmvnorm(g = g_star, Omega = params$Omega) + log(a_star) * p
  prop <- 
    exp(loglik_star - loglik)
  if(runif(1) < prop)
    return(list(a = a_star,
                loglik = loglik_star,
                accept = TRUE))
  else
    return(list(a = a,
                loglik = loglik,
                accept = FALSE))
}

a_to_u <- function(a, pi0, mu, sigma) {
  to_return <-  pi0 / 2
  to_return[a > 0] <- 
    pi0[a > 0] + 
    pnorm(log(a[a > 0]), 
          mean = mu[a > 0], 
          sd = sigma[a > 0]) * 
    (1 - pi0[a > 0])
  
  return(to_return)
}

log_dmvnorm <- function(g, Omega) {
  return((-t(g) %*% Omega %*% g / 2)[1, 1, drop = TRUE])
}

EM_fit <- function(data, initials, controls) {
  while(TRUE) {
    ## E step
    doParallel::registerDoParallel(cores = controls$ncores)
    eas <- foreach::`%dopar%`(
      foreach::foreach(i_sample = seq_len(nrow(data))),
      {
        a_samples <- SparseDOSSA2:::mcmc_a(x = data[i_sample, , drop = TRUE],
                                           params = params,
                                           R = controls$R)
        SparseDOSSA2:::one_step_a(x = data[i_sample, , drop = TRUE],
                                  params = params, a = a_samples[[80]]$a,
                                  loglik = a_samples[[80]]$a)
        return(mean(a_samples))
      })
    doParallel::stopImplicitCluster()
  }
}
