fit_featureParam <- function(data) {
  feature_param <- 
    t(vapply(seq_len(ncol(data)),
             function(i_feature)
               fit_oneFeatureParam(
                 data[, i_feature, drop = TRUE]),
             rep(0.0, 6)))
  # fill in missing sigmas (those with only one non-zero observations)
  # with the overall median across features
  feature_param[is.na(feature_param[, "sigma"]), "sigma"] <- 
    median(feature_param[, "sigma"], na.rm = TRUE)
  rownames(feature_param) <- colnames(data)
  return(feature_param)
}

fit_oneFeatureParam <- function(x) {
  ind_nonzero <- x > 0
  if(!any(ind_nonzero))
    stop("Does not support features with all zero values!")
  if(all(ind_nonzero)) ## FIXME
    pi0 <- 0.5 / length(x)
  else
    pi0 <- mean(!ind_nonzero)
  mu <- mean(log(x[ind_nonzero]))
  sigma <- sd(log(x[ind_nonzero]))
  
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
           "sigmaMod" = sigmaMod))
}

fit_F <- function(feature_param) {
  if(!all(colnames(feature_param) == c("pi0", "mu", "sigma")))
    stop("feature_param is not of the correct format!")
  # transform pi0 parameter before estimation
  feature_param[, "pi0"] <- log(feature_param[, "pi0"]) - log(1 - feature_param[, "pi0"])
  
  ind_zero_sigma <- feature_param[, "sigma"] == 0
  
  K_nonzero <- ks::Hscv(x = feature_param[!ind_zero_sigma, , drop = FALSE])
  if(any(ind_zero_sigma))
    ks::Hscv(x = feature_param[ind_zero_sigma, -3, drop = FALSE])
  else
    K_zero <- NULL
  
  return(list(p0_sigma = mean(ind_zero_sigma),
              K_nonzero = K_nonzero,
              K_zero = K_zero))
}

fit_depth <- function(depth) {
  return(c("mu_depth" = mean(log(depth)),
           "sigma_depth" = sd(log(depth))))
}
