#' Initial per-feature logit normal mean estimation
#'
#' @param y_ij vector of per-feature read count across samples
#' @param n_i vector of total read count across samples
#'
#' @return estimated log normal mean
get_muhat_i <- function(y_ij, n_i) {
  y_ij <- add_ifzero(y_ij, .5)
  phat_i <- y_ij/n_i
  log(phat_i) - log(1 - phat_i)
}

#' Initial per-feature logit normal variance estimation
#' @param y_ij vector of per-feature read count across samples
#' @param n_i vector of total read count across samples
#'
#' @return estimated log normal variance
get_sigma2hat_i <- function(y_ij, n_i) {
  y_ij <- add_ifzero(y_ij, 0.5)
  phat_i <- y_ij/n_i
  1 / (n_i * phat_i * (1 - phat_i))
}

#' Add a fixed value to zero counts (mainly used for Laplacian approximation of binomial likelihood)
#'
#' @param y_ij vector of per-feature read count across samples
#'
#' @return modified counts
add_ifzero <- function(y_ij, value) {
  ind_zero_i <- y_ij == 0
  y_ij[ind_zero_i] <- y_ij[ind_zero_i] + value
  y_ij
}


#' E step update for the per sample weights of prevalence parameter
#'
#' @param muhat_i vector of muhat estimation across samples
#' @param sigma2hat_i vector of sigma2hat estimation across samples
#' @param y_ij vector of per-feature read count across samples
#' @param n_i vector of total read count across samples
#' @param pi pi parameter (in current iteration)
#' @param mu mu parameter (in current iteration)
#' @param sigma2 sigma2 parameter (in current iteration)
#'
#' @return weights of observations, in the form as a) expectation (probability) of sequencing (as
#' opposed to biological) zero for zero-count samples, and b) ones for non-zero-count samples
get_w_i <- function(muhat_i, sigma2hat_i,
                    y_ij, n_i,
                    pi, mu, sigma2) {
  ind_zero_i <- y_ij == 0
  w_i <- rep(1, length = length(y_ij))
  likelihood_normal_i <- dnorm(muhat_i[ind_zero_i],
                               mean = mu,
                               sd = sqrt(sigma2hat_i[ind_zero_i] + sigma2))
  w_i[ind_zero_i] <- (pi * likelihood_normal_i) / (1 - pi + pi * likelihood_normal_i)
  w_i
}

#' M step update for pi
#'
#' @param w_i vector of observation weights
#'
#' @return updated pi estimation
update_pi <- function(w_i) {
  mean(w_i)
}

#' M step update for mu
#'
#' @param sigma2 sigma2 parameters
#' @param muhat_i vector of muhat estimation across samples
#' @param sigma2hat_i vector of sigma2hat estimation across samples
#' @param w_i vector of observation weights
#'
#' @return updated mu estimation
update_mu <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  weights_i <- w_i / (sigma2 + sigma2hat_i)
  sum(weights_i * muhat_i) / sum(weights_i)
}

#' M step update for sigma2
#'
#' @param muhat_i vector of muhat estimation across samples
#' @param sigma2hat_i vector of sigma2hat estimation across samples
#' @param w_i vector of observation weights
#' @param sigma2_old old sigma2 parameter (used as initial value for optimization)
#' @param control list of control parameters for the inner layer optim function
#'
#' @return updated sigma2 estimation
update_sigma2 <- function(muhat_i, sigma2hat_i, w_i, sigma2_old, control) {
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
        lower = 0, upper = Inf,
        control = control)$par
}

#' Negative log likelihood with respect to sigma2
#'
#' @param sigma2 sigma2 parameter
#' @param muhat_i vector of muhat estimation across samples
#' @param sigma2hat_i vector of sigma2hat estimation across samples
#' @param w_i vector of observation weights
#'
#' @return negative log likelihood, evaluated at sigma2
negLogLik_sigma2 <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  mu_sigma2_i <- update_mu(sigma2 = sigma2,
                           muhat_i = muhat_i,
                           sigma2hat_i = sigma2hat_i,
                           w_i = w_i)
  sum(w_i * log(sigma2 + sigma2hat_i)) +
    sum(w_i * (muhat_i - mu_sigma2_i)^2 / (sigma2 + sigma2hat_i))
}

#' Derivative of the negative log likelihood (for gradient decent)
#'
#' @param sigma2 sigma2 parameter
#' @param muhat_i vector of muhat estimation across samples
#' @param sigma2hat_i vector of sigma2hat estimation across samples
#' @param w_i vector of observation weights
#'
#' @return  derivative of negative log likelihood, evaluated at sigma2
dNegLogLik_sigma2 <- function(sigma2, muhat_i, sigma2hat_i, w_i) {
  mu_sigma2_i <- update_mu(sigma2 = sigma2,
                           muhat_i = muhat_i,
                           sigma2hat_i = sigma2hat_i,
                           w_i = w_i)
  sum(w_i / (sigma2 + sigma2hat_i)) -
    sum(w_i * (muhat_i - mu_sigma2_i)^2 / (sigma2 + sigma2hat_i)^2)

}

#' Posterior estimation of per-sample mu (shrinked by the estimated parameters)
#'
#' @param muhat_i vector of muhat estimation across samples (before shrinkage)
#' @param sigma2hat_i vector of sigma2hat estimation across samples (before shrinkage)
#' @param mu estimated mu parameter
#' @param sigma2 estimated sigma2 parameter
#'
#' @return posterior muhat_i, shrinked towards mu
posterior_mu_i <- function(muhat_i, sigma2hat_i, mu, sigma2) {
  (muhat_i * sigma2 + mu * sigma2hat_i) / (sigma2hat_i + sigma2)
}

#' Posterior estimation of per-sample sigma2 (shrinked by the estimated parameters)
#'
#' @param sigma2hat_i vector of sigma2hat estimation across samples (before shrinkage)
#' @param sigma2 estimated sigma2 parameter
#'
#' @return  posterior sigma2hat_i
posterior_sigma2_i <- function(sigma2hat_i, sigma2) {
  sigma2hat_i * sigma2 / (sigma2hat_i + sigma2)
}
