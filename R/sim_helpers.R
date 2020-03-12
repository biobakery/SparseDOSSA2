filter_params <- function(params, max_pi0 = 0.95) {
  ind <- params$pi0 < max_pi0
  return(list(pi0 = params$pi0[ind],
              mu = params$mu[ind],
              sigma = params$sigma[ind],
              Omega = threshold_matrix(solve(params$Sigma[ind, ind])),
              Sigma = params$Sigma[ind, ind]))
}

simulate_a <- function(n, params, k_feature = 2, maxit = 10, verbose = FALSE) {
  i_iter <- 0
  while(TRUE) {
    if(i_iter + 1 > maxit)
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

params_a_to_x <- function(params_a, n = 100000) {
  samples_a <- rcopulasso(n = n,
                          pi0 = params_a$pi0,
                          mu = params_a$mu,
                          sigma = params_a$sigma,
                          Omega = params_a$Omega)
  samples_x <- t(apply(samples_a, 1, function(x) x / sum(x)))
  mat_marginals <- get_marginals(samples_x)
  return(list(pi0 = mat_marginals[, "pi0"],
              mu = mat_marginals[, "mu"],
              sigma = mat_marginals[, "sigma"],
              Corr = cor(samples_x, method = "spear")))
}
