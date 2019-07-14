#' @param n_ij feature count of sample i feature j
#' @param n_i total read count of sample i
#' 
#' @return
#'
#' @examples
get_muhat_i <- function(n_ij, n_i) {
  n_ij <- addone_ifzero(n_ij)
  n_ij <- minusone_ifmax(n_ij, n_i) ## FIXME
  phat_i <- n_ij/n_i
  log(phat_i) - log(1 - phat_i)
}

#' @param n_ij feature count of sample i feature j
#' @param n_i total read count of sample i
#' 
#' @return
#'
#' @examples
get_sigma2hat_i <- function(n_ij, n_i) {
  n_ij <- addone_ifzero(n_ij)
  n_ij <- minusone_ifmax(n_ij, n_i) ## FIXME
  phat_i <- n_ij/n_i
  1 / (n_i * phat_i * (1 - phat_i))
}

#' @param n_ij feature count of sample i feature j
#' 
#' @return
#'
#' @examples
addone_ifzero <- function(n_ij) {
  ind_zero_i <- n_ij == 0
  n_ij[ind_zero_i] <- n_ij[ind_zero_i] + 1
  n_ij
}

#' @param n_ij feature count of sample i feature j
#' @param n_i total read count of sample i
#' 
#' @return
#'
#' @examples
minusone_ifmax <- function(n_ij, n_i) {
  ind_max_i <- n_ij >= n_i
  n_ij[ind_max_i] <- n_ij[ind_max_i] - 1
  n_ij
}

#' @param muhat_i mu iteration
#' @param sigma2hat_i sigma2 iteration
#' @param n_ij feature count of sample i feature j
#' @param n_i total read count of sample i
#' @param pi0 
#' @param mu 
#' @param sigma2 
#' 
#' @return weights
#'
#' @examples
get_w_i <- function(muhat_i, sigma2hat_i, 
                    n_ij, n_i,
                    pi0, mu, sigma2) {
  ind_zero_i <- n_ij == 0
  w_i <- rep(1, length = length(n_ij))
  likelihood_normal_i <- dnorm(muhat_i[ind_zero_i], 
                               mean = mu, 
                               sd = sqrt(sigma2hat_i[ind_zero_i] + sigma2))
  w_i[ind_zero_i] <- 1 - pi0 / (pi0 + (1 - pi0) * likelihood_normal_i)
  w_i
}

update_pi0 <- function(w_i) mean(1 - w_i)

update_mu <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  weights_i <- w_i / (sigma2 + sigma2hat_i)
  sum(weights_i * muhat_i) / sum(weights_i)
}

update_sigma2 <- function(muhat_i, sigma2hat_i, w_i, sigma2_old) {
  mu_sigma2is0_i <- update_mu(sigma2 = 0, 
                              muhat_i = muhat_i,
                              sigma2hat_i = sigma2hat_i,
                              w_i = w_i)
  if(sum(sum(w_i / sigma2hat_i^2 * (mu_sigma2is0_i - muhat_i)^2)) <= 
     sum(w_i / sigma2hat_i))
    return(0)
  optim(par = sigma2_old,
        fn = negLogLik_sigma2, 
        gr = dNegLogLik_sigma2,
        muhat_i = muhat_i,
        sigma2hat_i = sigma2hat_i,
        w_i = w_i,
        method = "L-BFGS-B",
        lower = 0, upper = Inf)$par
}

negLogLik_sigma2 <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  mu_sigma2_i <- update_mu(sigma2 = sigma2,
                           muhat_i = muhat_i,
                           sigma2hat_i = sigma2hat_i,
                           w_i = w_i)
  sum(w_i * log(sigma2 + sigma2hat_i)) + 
    sum(w_i * (muhat_i - mu_sigma2_i)^2 / (sigma2 + sigma2hat_i))
}

dNegLogLik_sigma2 <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  mu_sigma2_i <- update_mu(sigma2 = sigma2,
                           muhat_i = muhat_i,
                           sigma2hat_i = sigma2hat_i,
                           w_i = w_i)
  sum(w_i / (sigma2 + sigma2hat_i)) -
    sum(w_i * (muhat_i - mu_sigma2_i)^2 / (sigma2 + sigma2hat_i)^2)
  
}

posterior_mu_i <- function(muhat_i, sigma2hat_i, mu, sigma2) {
  (muhat_i * sigma2 + mu * sigma2hat_i) / (sigma2hat_i + sigma2) 
}

posterior_sigma2_i <- function(sigma2hat_i, sigma2) {
  sigma2hat_i * sigma2 / (sigma2hat_i + sigma2)
}
