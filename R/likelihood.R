pcopulasso <- function(x, 
                       mean, sd, pi0, sigma,
                       R = 100000) {
  samples_a <- SparseDOSSA2:::rcopulasso(n = R,
                                         mean = params$mu,
                                         sd = params$sigma,
                                         pi0 = params$pi0,
                                         sigma = solve(params$Omega))
  samples_asum <- vapply(seq_len(R), function(r) sum(samples_a[r, ]), 
                         0.0)
  
  for(i in seq_len(nrow(x))) {
    i_rx <- rx(c = x[i, ], C = sum(x[i, ]), R = R)
    samples_x <- i_rx$samples
    samples_xsum <- vapply(seq_len(R), function(r) sum(samples_a[r, ]), 
                           0.0)
  }
  l_rx <- lapply(seq_len(nrow(x)), 
                 function(i) {
                   rx(c = x[i, ], C = sum(x[i, ]), R = R)
                 })
  
  sum(log(vapply(seq_len(nrow(x)), 
                 function(i)
                   mean(vapply(seq_along(R), function(r) {
                     dmultinom(x = x[i, ], prob = samples_X[r, ])
                   }, 0.0)),
                 0.0)))
}

rx <- function(c, C, R) {
  x <- (c + 0.5) / (C + 0.5)
  mu <- log(x)
  sigma <- sqrt((1 - x) / (c + 0.5))
  
  samples <- exp(rnorm(n = R * length(c), mean = mu, sd = sigma)) * 
    rbinom(n = R * length(c), size = 1, prob = 0.5)
  samples <- t(matrix(samples, nrow = length(c)))
  params <- list(mu = mu, sigma = sigma)
  
  return(list(samples = samples, params = params))
}

a_to_u <- function(a, pi0, mu, sigma) {
  to_return <-  pi0 / 2
  if(any(a > 0)) ##FIXME
    to_return[a > 0] <- 
      pi0[a > 0] + 
      pnorm(log(a[a > 0]), 
            mean = mu[a > 0], 
            sd = sigma[a > 0]) * 
      (1 - pi0[a > 0])
  
  return(to_return)
}

dcopulasso <- function(x, mu, sigma, pi0, Omega, R)

u_to_g <- function(u, a, mu, sigma) {
  g <- qnorm(u)
  if(any(g == Inf)) ##FIXME
    g[g == Inf] <- (log(a[g == Inf]) - mu[g == Inf]) / sigma[g == Inf]
  
  return(g)
}

logLik_copula <- function(g, asum, x,
                          mu, sigma, 
                          Omega) {
  log_dmvnorm(S = g %*% t(g), Omega = Omega) + sum(g^2/2) - 
    sum((log(asum * x[x > 0]) - mu[x > 0])^2 / (sigma[x > 0])^2 / 2)
}

dmultinom2 <- function(x, prob) {
  
}