filter_params <- function(params, max_pi0 = 0.95) {
  ind <- params$pi0 < max_pi0
  Omega <- solve(params$Sigma[ind, ind])
  return(list(pi0 = params$pi0[ind],
              mu = params$pi0[ind],
              sigma = params$pi0[ind],
              Omega = solve(params$Sigma[ind, ind]),
              Sigma = params$Sigma[ind, ind]))
}

simulate_a <- function(n, params, k_feature = 2, max_iter = 10, verbose = FALSE) {
  i_iter <- 0
  while(TRUE) {
    if(i_iter + 1 > max_iter)
      stop("Can't satisfy conditions!")
    i_iter <- i_iter + 1
    
    samples_a <- rcopulasso(n = n,
                            pi0 = params$pi0,
                            mu = params$mu,
                            sigma = params$sigma,
                            Omega = params$Omega)
    
    if(all(apply(samples_a > 0, 2, sum) > 1)) {
      if(verbose) print(i_iter)
      return(samples_a)
    }
  }
}
