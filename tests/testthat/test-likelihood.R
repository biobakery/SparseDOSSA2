test_that("dloga works", {
  logas <- seq(-5, 5, by = 5e-3)
  gs <- qnorm(0.5 + 0.5 * pnorm(logas))
  coef_lm <- lm(gs[-(1:1200)] ~ logas[-(1:1200)])$coef
  
  dloga2 <- function(a,
                    pi0, mu, sigma, Omega, Sigma, 
                    log.p = TRUE) {
    
    u <- a_to_u(a,
                pi0 = pi0, mu = mu, sigma = sigma)
    g <- qnorm(u)
    g[a > 0] <- log(a[a > 0]) * coef_lm[2] + coef_lm[1]
    
    ind_nonzero <- a > 0
    if(any(abs(g) == Inf)) {
      log_d <- -Inf
    } else if(all(ind_nonzero)) {
      log_d <- 
        mvtnorm::dmvnorm(x = g,
                         mean = rep(0, length(g)),
                         sigma = Sigma,
                         log = TRUE) +
        sum(g^2) / 2 -
        sum((log(a) - mu)^2 / (sigma)^2 / 2) - 
        sum(log(sigma)) + 
        sum(log(1 - pi0))
    } else if(!any(ind_nonzero)) {
      log_d <- 
        log(mvtnorm::pmvnorm(
          lower = -Inf, 
          upper = g, 
          mean = rep(0, length(g)),
          sigma = Sigma))
    } else {
      log_d <- 
        log(mvtnorm::pmvnorm(
          lower = -Inf, 
          upper = g[!ind_nonzero], 
          mean = (-solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE],
                         Omega[!ind_nonzero, ind_nonzero, drop = FALSE]) %*% 
                    g[ind_nonzero])[, 1],
          sigma = solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE]))) +
        mvtnorm::dmvnorm(x = g[ind_nonzero],
                         mean = rep(0, length = sum(ind_nonzero)),
                         sigma = Sigma[ind_nonzero, ind_nonzero, drop = FALSE],
                         log = TRUE) +
        sum(g[ind_nonzero]^2) / 2 -
        sum((log(a[ind_nonzero]) - mu[ind_nonzero])^2 / (sigma[ind_nonzero])^2 / 2) - 
        sum(log(sigma[ind_nonzero])) + 
        sum(log(1 - pi0[ind_nonzero]))
    }
    
    if(log.p)
      return(log_d)
    else
      return(exp(log_d))
  }
  
  Sigma <- matrix(0.3, nrow = 4, ncol = 4)
  diag(Sigma) <- rep(1, 4)
  Omega <- solve(Sigma)
  
  tmpfun1 <- function(a, pi0, mu, sigma, Omega, Sigma) {
    dloga(a, pi0, mu, sigma, Omega, Sigma, log.p = FALSE) /
      prod(a[a > 0])
  }
  
  
  val1 <- cubature::hcubature(tmpfun1,
                      lowerLimit = c(0.000001, 0.000001, 0.000001, 0.000001),
                      upperLimit = c(100, 100, 100, 100),
                      pi0 = c(0.5, 0.5, 0.5, 0.5),
                      mu = c(0, 0, 0, 0),
                      sigma = c(1, 1, 1, 1),
                      Omega = Omega,
                      Sigma = Sigma,
                      tol = 1e-2)
  
  tmpfun2 <- function(a, pi0, mu, sigma, Omega, Sigma) {
    dloga(c(0, a), pi0, mu, sigma, Omega, Sigma, log.p = FALSE) /
      prod(a[a > 0])
  }
  
  val2 <- cubature::hcubature(tmpfun2,
                      lowerLimit = c(0.000001, 0.000001, 0.000001),
                      upperLimit = c(100, 100, 100),
                      pi0 = c(0.5, 0.5, 0.5, 0.5),
                      mu = c(0, 0, 0, 0),
                      sigma = c(1, 1, 1, 1),
                      Omega = Omega,
                      Sigma = Sigma,
                      tol = 1e-2) # should be 0.5^4
  
  
  tmpfun3 <- function(a, pi0, mu, sigma, Omega, Sigma) {
    dloga(c(0, 0, a), pi0, mu, sigma, Omega, Sigma, log.p = FALSE) /
      prod(a[a > 0])
  }
  
  val3 <- cubature::hcubature(tmpfun3,
                      lowerLimit = c(0.0000001, 0.0000001),
                      upperLimit = c(100, 100),
                      pi0 = c(0.5, 0.5, 0.5, 0.5),
                      mu = c(0, 0, 0, 0),
                      sigma = c(1, 1, 1, 1),
                      Omega = Omega,
                      Sigma = Sigma,
                      tol = 1e-2) # should be 0.5^4
  
  
  tmpfun4 <- function(a, pi0, mu, sigma, Omega, Sigma) {
    dloga(c(0, 0, 0, a), pi0, mu, sigma, Omega, Sigma, log.p = FALSE) /
      prod(a[a > 0])
  }
  
  val4 <- cubature::hcubature(tmpfun4,
                      lowerLimit = c(0.000001),
                      upperLimit = c(100),
                      pi0 = c(0.5, 0.5, 0.5, 0.5),
                      mu = c(0, 0, 0, 0),
                      sigma = c(1, 1, 1, 1),
                      Omega = Omega,
                      Sigma = Sigma,
                      tol = 1e-2) # should be 0.5^4
  
  
  val5 <- dloga(c(0, 0, 0, 0), pi0 = c(0.5, 0.5, 0.5, 0.5),
        mu = c(0, 0, 0, 0),
        sigma = c(1, 1, 1, 1),
        Omega = Omega,
        Sigma = Sigma,
        log.p = FALSE) # should be 0.5^4
  
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0),
                     mean = c(0, 0),
                     sigma = diag(c(1, 1)),
                     log = TRUE) -
      dloga(a = c(1, 1),
            pi0 = c(0, 0),
            mu = c(0, 0),
            sigma = c(1, 1),
            Omega = diag(c(1, 1)),
            Sigma = diag(c(1, 1)))  
  ),
  1e-15)
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0), 
                     mean = c(0, 0),
                     sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                     log = TRUE) -
      dloga_old(a = c(1, 1),
                pi0 = c(0, 0),
                mu = c(0, 0), 
                sigma = c(1, 1),
                Omega = solve(matrix(c(1, 0.5, 0.5, 1), 2, 2)))
  ),
  1e-15)
})


test_that("dloga_old works", {
  expect_lt(abs(
    # FIXME
    0.5 * 0.5 * dnorm(0) -
      dloga(a = c(1, 0),
            pi0 = c(0.5, 0.5),
            mu = c(0, 0),
            sigma = c(1, 1),
            Omega = diag(c(1, 1)),
            Sigma = diag(c(1, 1)), log.p = FALSE)  
  ),
  1e-15)
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0), 
                     mean = c(0, 0),
                     sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                     log = TRUE) -
      dloga_old(a = c(1, 1),
                pi0 = c(0, 0),
                mu = c(0, 0), 
                sigma = c(1, 1),
                Omega = solve(matrix(c(1, 0.5, 0.5, 1), 2, 2)))
  ),
  1e-15)
})

test_that("integrand_dx works", {
  expect_lt(abs(mvtnorm::dmvnorm(x = c(0, 0),
                                 mean = c(0, 0),
                                 sigma = diag(c(1, 1))) -
                  integrand_dx(log_asum = log(2),
                               x = c(0.5, 0.5),
                               mu = c(0, 0),
                               sigma = c(1, 1),
                               pi0 = c(0, 0),
                               Omega = diag(c(1, 1)),
                               log_offset = 0)
  ),
  1e-15)
  expect_lt(abs(mvtnorm::dmvnorm(x = c(0, 0),
                                 mean = c(0, 0),
                                 sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)) -
                  integrand_dx(log_asum = log(4),
                               x = c(0.25, 0.25),
                               mu = c(0, 0),
                               sigma = c(1, 1),
                               pi0 = c(0, 0),
                               Omega = solve(matrix(c(1, 0.5, 0.5, 1), 2, 2)),
                               log_offset = 0)
  ),
  1e-15)
})

test_that("dx_works", {
  dx_1d <- Vectorize(
    function(x, pi0, mu, sigma, Omega)
      dx(c(x, 1 - x),
         pi0, mu, sigma, Omega,
         control = list(jacobian = TRUE)),
    vectorize.args = "x")
  expect_lt(abs(integrate(dx_1d,
                          subdivisions = 10000,
                          rel.tol = .Machine$double.eps^0.5,
                          lower = 0, upper = 1,
                          pi0 = c(0, 0), mu = c(0, 0), sigma = c(1, 1),
                          Omega = diag(c(1, 1)))$value -
                  1),
            .Machine$double.eps^0.5)
})