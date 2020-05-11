estimate_featureParam_new <- function(x) {
  ind_nonzero <- x > 0
  if(sum(ind_nonzero) <= 1) ## FIXME
    stop("Does not support features with only one non-zero value or less!")
  if(all(!ind_nonzero)) ## FIXME
    stop("Does not support features with all non-zero values!")
  
  pi0 <- mean(!ind_nonzero)
  mu <- mean(log(x[ind_nonzero]))
  if(sum(ind_nonzero) > 1)
    sigma <- sd(log(x[ind_nonzero]))
  else
    sigma <- 1 ##FIXME
  
  # Additional parameters to help with computation
  glim <- qnorm(pi0)
  g0 <- truncnorm::etruncnorm(b = glim)
  sigmaMod <- sqrt(pi0 * g0^2 + 
                     (1 - pi0) * 
                     (truncnorm::vtruncnorm(a = glim) + 
                        truncnorm::etruncnorm(a = glim)^2))
  
  return(c("pi0" = pi0,
           "mu" = mu,
           "sigma" = sigma,
           "glim" = glim,
           "g0" = g0,
           sigmaMod = sigmaMod))
}

get_marginals <- function(X) {
  as.data.frame(
    t(vapply(seq_len(ncol(X)),
             function(i_feature)
               estimate_featureParam_new(
                 X[, i_feature, drop = TRUE]),
             rep(0.0, 6))))
}
