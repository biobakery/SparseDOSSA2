test_that("get_intLimits works", {
  f_test <- function(x, lim) {
    ifelse(abs(x) <= lim,
           1,
           0)
  }
  
  expect_equal(get_intLimits(f_test, 
                             limit_max = 2^10, limit_min = 1, step_size = 2,
                             lim = 2^5),
               c(-2^6, 2^6))
  expect_equal(get_intLimits(f_test, 
                             limit_max = 2^10, limit_min = 1, step_size = 2,
                             lim = 1),
               c(-2, 2))
  expect_warning(get_intLimits(f_test, 
                               limit_max = 2^10, limit_min = 2, step_size = 2,
                               lim = 1))
  expect_error(get_intLimits(f_test, 
                             limit_max = 2^10, limit_min = 2, step_size = 2,
                             lim = 2^11))
})