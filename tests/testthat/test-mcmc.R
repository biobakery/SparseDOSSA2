test_that("loglik_copula is correct", {
  expect_equal(loglik_copula(g = c(0, 0), 
                asum = 1, x = c(1, 0),
                mu = c(0, 0), sigma = c(1, 1), 
                Omega = diag(c(1, 1))),
               0)
})
