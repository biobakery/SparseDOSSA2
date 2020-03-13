estimate_featureParam_new <- function(x) {
  ind_nonzero <- x > 0
  if(all(ind_nonzero)) ## FIXME
    pi0 <- length(ind_nonzero) / (length(ind_nonzero) + 0.5)
  else
    pi0 <- mean(!ind_nonzero)
  mu <- mean(log(x[ind_nonzero]) - log(1 - x[ind_nonzero])) ## FIXME
  if(sum(ind_nonzero) > 1)
    sigma <- sd(log(x[ind_nonzero]) - log(1 - x[ind_nonzero]))
  else
    sigma <- 1 ##FIXME
  
  return(c("pi0" = pi0,
           "mu" = mu,
           "sigma" = sigma))
}

get_marginals <- function(X) {
  t(vapply(seq_len(ncol(X)),
           function(i_feature)
             estimate_featureParam_new(
               X[, i_feature, drop = TRUE]),
           c(0.0, 0.0, 0.0)))
}