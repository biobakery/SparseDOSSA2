estimate_readCount <- function(n_i) {
  return(c("mu" = mean(log(n_i)),
           "sigma2" = var(log(n_i))))
}