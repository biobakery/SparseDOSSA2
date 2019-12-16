test_that("dloga works", {
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0),
                     mean = c(0, 0),
                     sigma = diag(c(1, 1)),
                     log = TRUE) -
      dloga(a = Rmpfr::mpfr(c(1, 1), precBits = 100),
            pi0 = c(0, 0),
            mu = c(0, 0),
            sigma = c(1, 1),
            Omega = diag(c(1, 1)))  
  ),
  1e-15)
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0), 
                     mean = c(0, 0),
                     sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                     log = TRUE) -
      dloga(a = Rmpfr::mpfr(c(1, 1), precBits = 100),
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
                  integrand_dx(log_asum = log(Rmpfr::mpfr(2, precBits = 100)),
                               x = c(0.5, 0.5),
                               mu = c(0, 0),
                               sigma = c(1, 1),
                               pi0 = c(0, 0),
                               Omega = diag(c(1, 1)))
  ),
  1e-15)
  expect_lt(abs(mvtnorm::dmvnorm(x = c(0, 0),
                                 mean = c(0, 0),
                                 sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)) -
                  vintegrand_dx(log_asum = log(Rmpfr::mpfr(4, precBits = 100)),
                                x = c(0.25, 0.25),
                                mu = c(0, 0),
                                sigma = c(1, 1),
                                pi0 = c(0, 0),
                                Omega = solve(matrix(c(1, 0.5, 0.5, 1), 2, 2)))
  ),
  1e-15)
})

test_that("dx_works", {
  dx_1d <- function(x, pi0, mu, sigma, Omega) {
    Rmpfr::sapplyMpfr(x,
                      function(x) {
                        dx(c(x, 1 - x), 
                           pi0, mu, sigma, Omega,
                           control = list(jacobian = TRUE))
                      })
  }
  
  expect_lt(abs(Rmpfr::integrateR(
    dx_1d,
    lower = 0, upper = 1, 
    pi0 = c(0, 0), mu = c(0, 0), sigma = c(1, 1),
    Omega = diag(c(1, 1)))$value - 
      1),
    1e-15)
})