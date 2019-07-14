#' Title
#'
#' @param n_ij feature count of sample i feature j
#' @param n_i total read count of sample i
#' @param control 
#'
#' @return
#'
#' @examples
EM_theta <- function(n_ij, n_i, 
                     control = list(maxiter = 1000,
                                    threshold = 1e-9)) {
  # report suspicious input
  if(all(n_ij == 0))
    stop("All feature counts are zeros!")
  if(length(n_ij) != length(n_i))
    stop("Length of n_ij and n_i must agree!")
  if(any(n_i < 2))
    stop("n_i cannot be less than 2!")
  
  # zero indicator, muhat_i, and sigma2hat_i
  ind_zero_i <- n_ij == 0
  muhat_i <- get_muhat_i(n_ij = n_ij, n_i = n_i)
  sigma2hat_i <- get_sigma2hat_i(n_ij = n_ij, n_i = n_i)
  
  # Initial values for EM
  mu_origin <- mean(muhat_i[!ind_zero_i])
  if(sum(!ind_zero_i) == 1) {
    sigma2_origin <- 0
  } else {
    sigma2_origin <- var(muhat_i[!ind_zero_i])
  }
  pi0_origin <- mean(ind_zero_i)
  # mu_old <- 0
  # sigma2_old <- 1
  # pi0_old <- 0.5
  mu_old <- mu_origin
  sigma2_old <- sigma2_origin
  pi0_old <- pi0_origin
  counter <- 0
  
  # EM updates
  while(TRUE) {
    counter <- counter + 1
    if(counter > control$maxiter) {
      warning("Maximum iteration reached!")
      break
    }
    
    # E step
    w_i <- get_w_i(muhat_i = muhat_i, 
                   sigma2hat_i = sigma2hat_i,
                   n_ij = n_ij,
                   n_i = n_i,
                   pi0 = pi0_old,
                   mu = mu_old,
                   sigma2 = sigma2_old)
    
    # M step
    pi0_new <- update_pi0(w_i = w_i)
    sigma2_new <- update_sigma2(muhat_i = muhat_i,
                                sigma2hat_i = sigma2hat_i,
                                w_i = w_i,
                                sigma2_old = sigma2_old)
    mu_new <- update_mu(sigma2 = sigma2_new,
                        muhat_i = muhat_i,
                        sigma2hat_i = sigma2hat_i,
                        w_i = w_i)
    
    if(all(abs(c(mu_new, sigma2_new, pi0_new) -
               c(mu_old, sigma2_old, pi0_old)) < control$threshold)) break
    # update maximization iterations
    mu_old <- mu_new
    sigma2_old <- sigma2_new
    pi0_old <- pi0_new
  }
  
  mu_posterior_i <- posterior_mu_i(muhat_i, sigma2hat_i, 
                                   mu_new, sigma2_new)
  sigma2_posterior_i <- posterior_sigma2_i(sigma2hat_i, sigma2_new)
  
  return(list(theta = c(mu = mu_new,
                        sigma2 = sigma2_new,
                        pi0 = pi0_new),
              theta_original = c(mu = mu_origin,
                                 sigma2 = sigma2_origin,
                                 pi0 = pi0_origin),
              hidden_param = list(w_i = w_i,
                                  mu_sample_i = muhat_i,
                                  sigma2_sample_i = sigma2hat_i,
                                  mu_posterior_i = mu_posterior_i,
                                  sigma2_posterior_i = sigma2_posterior_i),
              niter = counter))
}
