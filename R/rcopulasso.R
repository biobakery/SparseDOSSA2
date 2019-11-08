rcopulasso <- function(n, mean, sd, pi0, sigma) {
  if(length(mean) != length(sd) |
     length(mean) != length(pi0) |
     length(mean) != nrow(sigma))
    stop("Parameter dimensions must agree!")
  
  # sample marginals
  mat_amarginals <- 
    vapply(seq_len(length(mean)),
           function(i)
             rZILogN_one(n = n,
                           mean = mean[i],
                           sd = sd[i],
                           pi0 = pi0[i]),
           rep(0.0, n))
  # arrange from smallest to largest for shuffling
  mat_amarginals <- 
    apply(mat_amarginals, 2, function(x) x[order(x)])
  
  # sample ranks
  mat_rank <- 
    mvtnorm::rmvnorm(n = n, sigma = sigma)
  mat_rank <- apply(mat_rank, 2, rank)
  
  mat_a <- 
    vapply(seq_len(length(mean)),
           function(i)
             mat_amarginals[, i, drop = TRUE][mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  
  return(mat_a)
}

rZILogN_one <- function(n, mean, sd, pi0) {
  return(exp(rnorm(n = n, mean = mean, sd = sd)) *
           rbinom(n = n, size = 1, prob = 1 - pi0))
}