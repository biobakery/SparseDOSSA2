estimate_featureParam_new <- function(x) {
  ind_nonzero <- x > 0
  pi0 <- mean(!ind_zero)
  mu <- mean(log(x[ind_nonzero]))
  if(sum(ind_nonzero) > 1)
    sigma <- sd(sd(log(x[ind_nonzero])))
  else
    sigma <- 1 ##FIXME

  return(c("pi0" = pi0,
           "mu" = mu,
           "sigma" = sigma))
}
