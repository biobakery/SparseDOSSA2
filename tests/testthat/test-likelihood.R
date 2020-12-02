test_that("dloga works", {
  pi0 <- rep(0.5, 3)
  mu <- rep(0, 3)
  sigma <- rep(1, 3)
  Sigma <- diag(rep(1, 3))
  Omega <- solve(Sigma)
  
  # certain dloga values should be the same as directly caculated from dnorm
  expect_equal(dnorm(0)^2 * 0.125,
               dloga(a = c(0, 1, 1),
                     pi0 = pi0,
                     mu = mu, 
                     sigma = sigma,
                     Omega = Omega, 
                     Sigma = Sigma, 
                     log.p = FALSE))
  
  int_func <- function(x) dloga(a = c(0, exp(x)),
                                pi0 = pi0,
                                mu = mu, 
                                sigma = sigma,
                                Omega = Omega, 
                                Sigma = Sigma, 
                                log.p = FALSE)
  
  # dloga should integrate out to discrete components of the likelihood
  expect_lt(abs(0.125 - 
                  cubature::pcubature(int_func, lowerLimit = rep(-10, 2), upperLimit = rep(10, 2))$integral),
            1e-5)
  
  # test for non trivial correlation
  Sigma <- matrix(0.1, 3, 3)
  diag(Sigma) <- rep(1, 3)
  Omega <- solve(Sigma)
  
  int_func <- function(x) dloga(a = c(0, exp(x)),
                                pi0 = pi0,
                                mu = mu, 
                                sigma = sigma,
                                Omega = Omega, 
                                Sigma = Sigma, 
                                log.p = FALSE)
  
  expect_lt(
    abs(
    cubature::pcubature(f = int_func, lower = rep(-10, 2), upper = rep(10, 2))$integral -
    mvtnorm::pmvnorm(lower = c(-Inf, 0, 0), upper = c(0, Inf, Inf), 
                     mean = rep(0, 3), sigma = Sigma)),
    1e-5)
})