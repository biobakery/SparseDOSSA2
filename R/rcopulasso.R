rcopulasso <- function(n, pi0, mu, sigma, Omega) {
  if(length(pi0) != length(mu) |
     length(pi0) != length(sigma) |
     length(pi0) != nrow(Omega))
    stop("Parameter dimensions must agree!")
  
  # sample marginals
  mat_amarginals <- 
    vapply(seq_len(length(pi0)),
           function(i)
             rZILogN_one(n = n,
                         pi0 = pi0[i],
                         mu = mu[i],
                         sigma = sigma[i]),
           rep(0.0, n))
  # arrange from smallest to largest for shuffling
  mat_amarginals <- 
    apply(mat_amarginals, 2, function(x) x[order(x)])
  
  # sample ranks
  mat_rank <- 
    mvtnorm::rmvnorm(n = n, sigma = solve(Omega))
  mat_rank <- apply(mat_rank, 2, rank)
  
  mat_a <- 
    vapply(seq_len(length(pi0)),
           function(i)
             mat_amarginals[, i, drop = TRUE][mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  
  return(mat_a)
}

rZILogN_one <- function(n, pi0, mu, sigma) {
  return(exp(rnorm(n = n, mean = mu, sd = sigma)) *
           rbinom(n = n, size = 1, prob = 1 - pi0))
}