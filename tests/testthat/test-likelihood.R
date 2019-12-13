test_that("dloga works", {
  expect_lt(abs(
    mvtnorm::dmvnorm(x = c(0, 0),
                     mean = c(0, 0),
                     sigma = diag(c(1, 1)),
                     log = TRUE) -
      dloga(a = c(1, 1),
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
      dloga(a = c(1, 1),
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
                               Omega = diag(c(1, 1)))
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
                               Omega = solve(matrix(c(1, 0.5, 0.5, 1), 2, 2)))
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