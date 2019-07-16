#' EM algorithm estimator for per-feature parameters (prevalence, mean, variance)
#'
#' @param y_ij vector of per-feature read count across samples
#' @param n_i vector of total read count across samples
#' @param control list of control parameters of EM algorithm
#'
#' @return list with components for estimated parameters, hidden parameters, and optimization information
estimate_featureParam <- function(y_ij, n_i,
                                  control = list(maxiter.outer = 100,
                                                 maxiter.inner = 1000,
                                                 reltol.outer = 1e-8,
                                                 factr.inner = 1e7)) {
  # sanity check on input
  if(all(y_ij == 0))
    stop("All feature counts are zeros!")
  if(length(y_ij) != length(n_i))
    stop("Length of y_ij and n_i must agree!")
  if(any(n_i < 2))
    stop("n_i cannot be less than 2!")

  # non-zero indicator, muhat_ij, and sigma2hat_ij
  ind_i <- y_ij != 0
  muhat_ij <- get_muhat_ij(y_ij = y_ij, n_i = n_i)
  sigma2hat_ij <- get_sigma2hat_ij(y_ij = y_ij, n_i = n_i)
  muhat_ij[!ind_i] <- muhat_ij[!ind_i] - sigma2hat_ij[!ind_i]

  # Initial values for EM
  mu_origin <- mean(muhat_ij[ind_i])
  if(sum(ind_i) == 1) {
    sigma2_origin <- 0
  } else {
    sigma2_origin <- var(muhat_ij[ind_i])
  }
  pi_origin <- mean(ind_i)

  mu_old <- mu_origin
  sigma2_old <- sigma2_origin
  pi_old <- pi_origin
  counter <- 0

  # EM updates
  while(TRUE) {
    counter <- counter + 1
    if(counter > control$maxiter.outer) {
      warning("Maximum iteration reached!")
      break
    }

    # E step
    w_i <- get_w_i(muhat_ij = muhat_ij,
                   sigma2hat_ij = sigma2hat_ij,
                   y_ij = y_ij,
                   n_i = n_i,
                   pi = pi_old,
                   mu = mu_old,
                   sigma2 = sigma2_old)

    # M step
    pi_new <- update_pi(w_i = w_i)
    sigma2_new <- update_sigma2(muhat_ij = muhat_ij,
                                sigma2hat_ij = sigma2hat_ij,
                                w_i = w_i,
                                sigma2_old = sigma2_old,
                                control = list(maxit = control$maxiter.inner,
                                               factr = control$factr.inner))
    mu_new <- update_mu(sigma2 = sigma2_new,
                        muhat_ij = muhat_ij,
                        sigma2hat_ij = sigma2hat_ij,
                        w_i = w_i)

    if(sigma2_old != 0) {
      if(all(abs((c(mu_new, sigma2_new, pi_new) -
                  c(mu_old, sigma2_old, pi_old)) /
                 c(mu_old, sigma2_old, pi_old)) <
             control$reltol.outer)) break
    } else {
      if(all(c(abs((c(mu_new, pi_new) -
                    c(mu_old, pi_old)) /
                   c(mu_old, pi_old)), sigma2_new) <
             control$reltol.outer)) break
    }

    # update maximization iterations
    mu_old <- mu_new
    sigma2_old <- sigma2_new
    pi_old <- pi_new
  }

  mu_posterior_ij <- posterior_mu_ij(muhat_ij, sigma2hat_ij,
                                     mu_new, sigma2_new)
  sigma2_posterior_ij <- posterior_sigma2_ij(sigma2hat_ij, sigma2_new)

  return(list(theta = c(mu = mu_new,
                        sigma2 = sigma2_new,
                        pi = pi_new),
              theta_original = c(mu = mu_origin,
                                 sigma2 = sigma2_origin,
                                 pi = pi_origin),
              hidden_param = list(w_i = w_i,
                                  mu_original_ij = muhat_ij,
                                  sigma2_original_ij = sigma2hat_ij,
                                  mu_posterior_ij = mu_posterior_ij,
                                  sigma2_posterior_ij = sigma2_posterior_ij),
              niter = counter))
}

#' Non-parametric density estimator for the distribution of per-feature parameters
#'
#' @param feature_params data frame of estimated per-feature parameters (in order mu, sigma2, and pi)
#' @param control list of control parameters to pass on to the bandwidth estimators in ks
#'
#' @return list with components for bandwidth estimations for zero and non-zero sigma2 part, and
#' percentage of zero sigma2s
#' @import ks
estimate_F <- function(feature_params, control = list(estimator = "Hscv")) {
  if(!all(colnames(feature_params) == c("mu", "sigma2", "pi")))
    stop("feature_params is not of the correct format!")
  # transform pi parameter before estimation
  feature_params[, 3] <- log(feature_params[, 3]) - log(1 - feature_params[, 3])

  ind_zero_sigma2 <- feature_params[, 2] == 0

  K_nonzero <- do.call(control$estimator, list(x = feature_params[!ind_zero_sigma2, , drop = FALSE]))
  if(any(ind_zero_sigma2))
    K_zero <- do.call(control$estimator, list(x = feature_params[ind_zero_sigma2, -2, drop = FALSE]))
  else
    K_zero <- NULL
  return(list(p0_sigma2 = mean(ind_zero_sigma2),
              K_nonzero = K_nonzero,
              K_zero = K_zero))
}

#' Normal copula estimator for feature-feature association (correlation between the p_ijs)
#'
#' @param feature_abd feature x sample abundance table
#' @param seed random seed (used for randomized correlation estimation)
#' @param control additional control parameters for correlation estimation
#'
#' @return the estimated copula
estimate_C <- function(feature_abd, seed,
                       control = list(method = "spearman",
                                      random = TRUE,
                                      R = 50)) {
  mat_cor <- do.call("cor2", c(list(x = t(feature_abd), seed = seed), control))
  copula::normalCopula(param = copula::P2p(mat_cor),
                       dim = nrow(feature_abd),
                       dispstr = "un")
}

estimate_readCount <- function(n_i) {
  return(c("mu" = mean(log(n_i)),
           "sigma2" = var(log(n_i))))
}

