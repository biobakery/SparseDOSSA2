#' EM algorithm estimator for per-feature parameters (prevalence, mean, variance)
#'
#' @param y_ij vector of per-feature read count across samples
#' @param n_i vector of total read count across samples
#' @param control list of control parameters of EM algorithm
#'
#' @return list with components for estimated parameters, hidden parameters, and optimization information
estimate_featureParams <- function(y_ij, n_i,
                                   control = list(maxiter.outer = 100,
                                                  maxiter.inner = 1000,
                                                  reltol.outer = 1e-8,
                                                  reltol.inner = 1e-5)) {
  # sanity check on input
  if(all(y_ij == 0))
    stop("All feature counts are zeros!")
  if(length(y_ij) != length(n_i))
    stop("Length of y_ij and n_i must agree!")
  if(any(n_i < 2))
    stop("n_i cannot be less than 2!")

  # non-zero indicator, muhat_ij, and sigma2hat_ij
  ind_i <- y_ij != 0
  muhat_ij <- get_muhat_i(y_ij = y_ij, n_i = n_i)
  sigma2hat_ij <- get_sigma2hat_i(y_ij = y_ij, n_i = n_i)

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
                                               reltol = control$reltol.inner))
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

  mu_posterior_i <- posterior_mu_i(muhat_ij, sigma2hat_ij,
                                   mu_new, sigma2_new)
  sigma2_posterior_i <- posterior_sigma2_i(sigma2hat_ij, sigma2_new)

  return(list(theta = c(mu = mu_new,
                        sigma2 = sigma2_new,
                        pi = pi_new),
              theta_original = c(mu = mu_origin,
                                 sigma2 = sigma2_origin,
                                 pi = pi_origin),
              hidden_param = list(w_i = w_i,
                                  mu_original_i = muhat_ij,
                                  sigma2_original_i = sigma2hat_ij,
                                  mu_posterior_i = mu_posterior_i,
                                  sigma2_posterior_i = sigma2_posterior_i),
              niter = counter))
}

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

estimate_readCount <- function(n_i) {
  return(c("mu" = mean(log(n_i)),
           "sigma2" = var(log(n_i))))
}

