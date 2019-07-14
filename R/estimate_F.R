estimate_F <- function(feature_params) {
  pi0_sigma2 <- mean(feature_params[, 2] == 0)
  K_nonzero <- ks::Hscv(x = feature_params[feature_params[, 2] > 0, ])
  if(pi0_sigma2 > 0)
    K_zero <- ks::Hscv(x = feature_params[feature_params[, 2] == 0, -2])
  else
    K_zero <- NULL
  return(list(pi0_sigma2 = pi0_sigma2, 
              K_nonzero = K_nonzero, 
              K_zero = K_zero))
}
