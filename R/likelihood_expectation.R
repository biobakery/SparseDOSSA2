dx <- function(x, 
               pi0, mu, sigma, Omega, Sigma,
               control = list(),
               log.p = FALSE) {
  control <- do.call(control_integrate, control)
  limits <- get_intLimits(
    x = x,
    pi0 = pi0, mu = mu, sigma = sigma, 
    Omega = Omega, Sigma = Sigma, 
    maxit = control$maxit_limits)
    
  fit_integrate <- 
    integrate2(vintegrand_dx,
               lower = limits[1], upper = limits[2], 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               precBits = control$precBits,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  # jacobian
  fit_integrate$integral <- as.double(log(fit_integrate$integral) - sum(log(x[x > 0])))
  fit_integrate$error <- as.double(exp(log(fit_integrate$error) -  sum(log(x[x > 0]))))
  
  if(log.p) {
    return(fit_integrate$integral)
  }
  else {
    fit_integrate$integral <- exp(fit_integrate$integral)
    if(control$only_value)
      return(fit_integrate$integral)
    return(fit_integrate)
  }
}

#' control parameters for the 
#' numerical integrations during the E step of SparseDOSSA2's fitting
#'
#' @param rel_tol relative change threshold in the integration values for the 
#' integration to converge
#' @param abs_tol absolute change threshold in the integration values for the 
#' integration to converge
#' @param max_eval maximum of integration evaluations allowed
#' @param maxit_limits maximum number of tries allowed to guess the integration's
#' lower and upper limits
#' @param precBits numeric precision used for the integration values
#' @param only_value whether or not only the integration value should be returned
#'
#' @return a list of the same names
#' @export
control_integrate <- function(rel_tol = 1e-2,
                              abs_tol = 0,
                              max_eval = 50,
                              maxit_limits = 10,
                              precBits = 200,
                              only_value = TRUE) {
  list(rel_tol = rel_tol,
       abs_tol = abs_tol,
       max_eval = max_eval,
       maxit_limits = maxit_limits,
       precBits = precBits,
       only_value = only_value)
}

ploga <- function(a, 
                  pi0, mu, sigma, Sigma,
                  log.p = TRUE) {
  ind_zero <- a == 0
  
  if(!any(ind_zero)) {
    log_p <- 1
  } else {
    a <- a[ind_zero]
    pi0 <- pi0[ind_zero]
    mu <- mu[ind_zero]
    sigma <- sigma[ind_zero]
    Sigma <- Sigma[ind_zero, ind_zero, drop = FALSE]
    
    u <- a_to_u(a,
                pi0 = pi0, mu = mu, sigma = sigma)
    g <- qnorm(u)
    log_p <- 
      log(mvtnorm::pmvnorm(
        lower = -Inf, 
        upper = g, 
        mean = rep(0, length(g)),
        sigma = Sigma))
  }
  
  if(log.p)
    return(log_p)
  else
    return(exp(log_p))
}

mloga <- function(a, 
                  pi0, mu, sigma, Sigma,
                  doComputeVariance = FALSE,
                  log_ploga = 0) {
  ind_zero <- a == 0
  
  if(!any(ind_zero)) {
    return(list(mean_cond = NULL,
                Sigma_cond = matrix(nrow = 0, ncol = 0)))
  } else {
    a <- a[ind_zero]
    pi0 <- pi0[ind_zero]
    mu <- mu[ind_zero]
    sigma <- sigma[ind_zero]
    Sigma <- Sigma[ind_zero, ind_zero, drop = FALSE]
    
    u <- a_to_u(a,
                pi0 = pi0, mu = mu, sigma = sigma)
    g <- qnorm(u)
    m_cond <- 
      tmvtnorm::mtmvnorm(sigma = Sigma,
                         upper = g,
                         doComputeVariance = doComputeVariance)
  }
  
  return(list(mean_cond = m_cond$tmean,
              Sigma_cond = m_cond$tvar))
}

dloga_forInt <- function(a,
                         pi0, mu, sigma, Omega, Sigma, 
                         log_ploga = 0,
                         mean_cond,
                         Sigma_cond = NA,
                         mean_cond_choice = "half",
                         Sigma_cond_choice = "cond",
                         log.p = TRUE) {
  
  u <- a_to_u(a,
              pi0 = pi0, mu = mu, sigma = sigma,
              half_pi0 = TRUE)
  g <- qnorm(u)
  
  ind_nonzero <- a > 0
  if(any(abs(g) == Inf)) {
    log_d <- -Inf
  } else if(all(ind_nonzero)) {
    log_d <- 
      mvtnorm::dmvnorm(x = g,
                       mean = rep(0, length(g)),
                       sigma = Sigma,
                       log = TRUE) +
      sum(g^2) / 2 -
      sum((log(a) - mu)^2 / (sigma)^2 / 2) - 
      sum(log(sigma)) + 
      sum(log(1 - pi0))
  } else if(!any(ind_nonzero)) {
    log_d <- 
      log_ploga
  } else {
    if(mean_cond_choice == "half")
      mean_zero <- g[!ind_nonzero]
    else
      mean_zero <- mean_cond
    if(Sigma_cond_choice == "cond")
      Sigma_nonzero <- solve(Omega[ind_nonzero, ind_nonzero, drop = FALSE])
    else if (Sigma_cond_choice == "full")
      Sigma_nonzero <- Sigma[ind_nonzero, ind_nonzero, drop = FALSE]
    else
      Sigma_nonzero <- solve(Omega[ind_nonzero, ind_nonzero, drop = FALSE]) +
        Sigma[ind_nonzero, !ind_nonzero, drop = FALSE] %*%
        solve(Sigma[!ind_nonzero, !ind_nonzero, drop = FALSE]) %*%
        Sigma_cond %*%
        solve(Sigma[!ind_nonzero, !ind_nonzero, drop = FALSE]) %*%
        Sigma[!ind_nonzero, ind_nonzero, drop = FALSE]
    log_d <- 
      log_ploga +
      mvtnorm::dmvnorm(
        x = g[ind_nonzero],
        mean = (Sigma[ind_nonzero, !ind_nonzero, drop = FALSE] %*%
                  solve(Sigma[!ind_nonzero, !ind_nonzero, drop = FALSE],
                        mean_zero))[, 1],
        sigma = Sigma_nonzero,
        log = TRUE) +
      sum(g[ind_nonzero]^2) / 2 -
      sum((log(a[ind_nonzero]) - mu[ind_nonzero])^2 / 
            (sigma[ind_nonzero])^2 / 2) - 
      sum(log(sigma[ind_nonzero])) + 
      sum(log(1 - pi0[ind_nonzero]))
  }
  
  if(log.p)
    return(log_d)
  else
    return(exp(log_d))
}

dloga <- function(a,
                  pi0, mu, sigma, Omega, Sigma, 
                  log.p = TRUE) {
  
  u <- a_to_u(a,
              pi0 = pi0, mu = mu, sigma = sigma)
  g <- qnorm(u)
  
  ind_nonzero <- a > 0
  if(any(abs(g) == Inf)) {
    log_d <- -Inf
  } else if(all(ind_nonzero)) {
    log_d <- 
      mvtnorm::dmvnorm(x = g,
                       mean = rep(0, length = length(pi0)),
                       sigma = Sigma,
                       log = TRUE) +
      sum(g^2) / 2 -
      sum((log(a) - mu)^2 / 
            sigma^2 / 2) - 
      sum(log(sigma)) + 
      sum(log(1 - pi0))
  } else if(!any(ind_nonzero)) {
    log_d <- 
      pmvnorm2(
        lower = -Inf, 
        upper = g, 
        mean = rep(0, length(g)),
        sigma = Sigma, 
        log.p = TRUE)
  } else {
    log_d <- 
      pmvnorm2(
        lower = -Inf, 
        upper = g[!ind_nonzero], 
        mean = (-solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE],
                       Omega[!ind_nonzero, ind_nonzero, drop = FALSE]) %*% 
                  g[ind_nonzero])[, 1],
        sigma = solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE]),
        log.p = TRUE) +
      mvtnorm::dmvnorm(x = g[ind_nonzero],
                       mean = rep(0, length = sum(ind_nonzero)),
                       sigma = Sigma[ind_nonzero, ind_nonzero, drop = FALSE],
                       log = TRUE) +
      sum(g[ind_nonzero]^2) / 2 -
      sum((log(a[ind_nonzero]) - mu[ind_nonzero])^2 / 
            (sigma[ind_nonzero])^2 / 2) - 
      sum(log(sigma[ind_nonzero])) + 
      sum(log(1 - pi0[ind_nonzero]))
  }
  
  if(log.p)
    return(log_d)
  else
    return(exp(log_d))
}

pmvnorm2 <- function(lower = -Inf, upper = Inf, 
                     mean, sigma, log.p = FALSE) {
  if(length(lower) == 1)
    lower <- rep(lower, times = length(mean))
  if(length(upper) == 1)
    upper <- rep(upper, times = length(mean))
  
  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = (abs(sigma) > 0) * 1,
    mode = "undirected",
    diag = FALSE
  )
  comp <- igraph::components(g)
  vp <- 
    vapply(seq_len(comp$no),
           function(i_comp) {
             i_ind <- comp$membership == i_comp
             log(mvtnorm::pmvnorm(lower = lower[i_ind],
                                  upper = upper[i_ind],
                                  mean = mean[i_ind],
                                  sigma = sigma[i_ind, i_ind, drop = FALSE]))
           },
           0.0)
  logp <- sum(vp)
  if(log.p)
    return(logp)
  else
    return(exp(logp))
}

dloga_old <- function(a,
                      pi0, mu, sigma, Omega,
                      log.p = TRUE) {
  u <- a_to_u_old(a,
                  pi0 = pi0, mu = mu, sigma = sigma)
  g <- qnorm(u)
  
  if(any(abs(g) == Inf)) {
    log_d <- -Inf
  } else {
    log_d <- du(g, Omega) -
      sum((log(a[a > 0]) - mu[a > 0])^2 / (sigma[a > 0])^2 / 2) -
      log(2 * pi) / 2 * sum(a > 0) - sum(log(sigma[a > 0]))
  }
  
  if(log.p)
    return(log_d)
  else
    return(exp(log_d))
}

du <- function(g, Omega, log.p = TRUE) {
  # Without normalizing constant (2pi)^(-2/p)!
  if(any(g == -Inf | g == Inf)) return(-Inf)
  log_d <- log_dmvnorm(S = g %*% t(g), Omega = Omega) + sum(g^2)/2 +
    log(det(Omega)) / 2
  
  if(log.p)
    return(log_d)
  else
    return(exp(log_d))
}

log_dmvnorm <- function(S, Omega) {
  - sum(Omega * S) / 2
}

integrand_dx <- function(log_asum, x, 
                         pi0, mu, sigma, 
                         Omega, Sigma) {
  dloga(a = a(x, exp(log_asum)),
        pi0 = pi0, mu = mu, sigma = sigma, 
        Omega = Omega, Sigma = Sigma,
        log.p = FALSE)
}

vintegrand_dx <- Vectorize2(integrand_dx, vectorize.args = "log_asum")

ea <- function(x, 
               pi0, mu, sigma, Omega, Sigma,
               control = list()) {
  control <- do.call(control_integrate, control)
  limits <- get_intLimits(
    x = x,
    pi0 = pi0, mu = mu, sigma = sigma, 
    Omega = Omega, Sigma = Sigma,
    maxit = control$maxit_limits)
  
  fit_integrate <- 
    integrate2(vintegrand_ea,
               lower = limits[1], upper = limits[2], 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               precBits = control$precBits,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  # jacobian
  fit_integrate$integral <- as.double(log(fit_integrate$integral) - sum(log(x[x > 0])))
  fit_integrate$error <- as.double(exp(log(fit_integrate$error) -  sum(log(x[x > 0]))))
  
  if(log.p) {
    return(fit_integrate$integral)
  }
  else {
    fit_integrate$integral <- exp(fit_integrate$integral)
    if(control$only_value)
      return(fit_integrate$integral)
    return(fit_integrate)
  }
}

integrand_ea <- function(log_asum, x, 
                         pi0, mu, sigma, 
                         Omega, Sigma) {
  exp(dloga(a = a(x, exp(log_asum)),
            pi0 = pi0, mu = mu, sigma = sigma, 
            Omega = Omega, Sigma = Sigma,
            log.p = TRUE) + 
        log_asum)
}

vintegrand_ea <- Vectorize2(integrand_ea, 
                            vectorize.args = "log_asum")

eloga <- function(x, 
                  pi0, mu, sigma, Omega, Sigma,
                  control = list()) {
  control <- do.call(control_integrate, control)
    limits <- get_intLimits(
      x = x,
      pi0 = pi0, mu = mu, sigma = sigma, 
      Omega = Omega, Sigma = Sigma,
      maxit = control$maxit_limits)
  
  fit_integrate <- 
    integrate2(vintegrand_eloga,
               lower = limits[1], upper = limits[2], 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               precBits = control$precBits,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  # jacobian
  fit_integrate$integral <- as.double(log(fit_integrate$integral) - sum(log(x[x > 0])))
  fit_integrate$error <- as.double(exp(log(fit_integrate$error) -  sum(log(x[x > 0]))))
  
  if(log.p) {
    return(fit_integrate$integral)
  }
  else {
    fit_integrate$integral <- exp(fit_integrate$integral)
    if(control$only_value)
      return(fit_integrate$integral)
    return(fit_integrate)
  }
}

integrand_eloga <- function(log_asum, x, 
                            pi0, mu, sigma, 
                            Omega, Sigma) {
  dloga(a = a(x, exp(log_asum)),
        pi0 = pi0, mu = mu, sigma = sigma, 
        Omega = Omega, Sigma = Sigma,
        log.p = FALSE) * log_asum
}

vintegrand_eloga <- Vectorize2(integrand_eloga, 
                               vectorize.args = "log_asum")

eloga2 <- function(x, 
                   pi0, mu, sigma, Omega, Sigma,
                   control = list()) {
  control <- do.call(control_integrate, control)
    limits <- get_intLimits(
      x = x,
      pi0 = pi0, mu = mu, sigma = sigma, 
      Omega = Omega, Sigma = Sigma,
      maxit = control$maxit_limits)
    
  fit_integrate <- 
    integrate2(vintegrand_eloga2,
               lower = limits[1], upper = limits[2], 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               precBits = control$precBits,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  # jacobian
  fit_integrate$integral <- as.double(log(fit_integrate$integral) - sum(log(x[x > 0])))
  fit_integrate$error <- as.double(exp(log(fit_integrate$error) -  sum(log(x[x > 0]))))
  
  if(log.p) {
    return(fit_integrate$integral)
  }
  else {
    fit_integrate$integral <- exp(fit_integrate$integral)
    if(control$only_value)
      return(fit_integrate$integral)
    return(fit_integrate)
  }
}

integrand_eloga2 <- function(log_asum, x, 
                             pi0, mu, sigma, 
                             Omega, Sigma) {
  dloga(a = a(x, exp(log_asum)),
        pi0 = pi0, mu = mu, sigma = sigma, 
        Omega = Omega, Sigma = Sigma,
        log.p = FALSE) * log_asum^2
}

vintegrand_eloga2 <- Vectorize2(integrand_eloga2, 
                                vectorize.args = "log_asum")

get_es <- function(x, pi0, mu, sigma, Omega, Sigma,
                   control) {
  time_start <- Sys.time()
  control <- do.call(control_integrate, control)
  
  limits <- get_intLimits(
    x = x,
    pi0 = pi0, mu = mu, sigma = sigma, 
    Omega = Omega, Sigma = Sigma,
    maxit = control$maxit_limits)
  
  neval <- 2
  knots_spline <- Rmpfr::mpfr(c(limits[1], limits[2]),
                              precBits = control$precBits)
  vals_spline <- 
    Rmpfr::mpfr(
      vintegrand_dx(as.double(knots_spline), 
                    pi0 = pi0, x = x, mu = mu, sigma = sigma,
                    Omega = Omega, Sigma = Sigma), 
      precBits = control$precBits)
  errors_spline <- Inf
  
  # find knots using dx
  while(TRUE) {
    i_max_error <- which(errors_spline == max(errors_spline))[1]
    knots_spline <- c(knots_spline[seq(1, i_max_error)],
                      Rmpfr::mean(knots_spline[c(i_max_error, i_max_error + 1)]),
                      knots_spline[seq(i_max_error + 1, neval)])
    vals_spline <- c(vals_spline[seq(1, i_max_error)],
                     Rmpfr::mpfr(vintegrand_dx(as.double(knots_spline[i_max_error + 1])),
                                 precBits = control$precBits),
                     vals_spline[seq(i_max_error + 1, neval)])
    
    neval <- neval + 1
    knots_diff <-  knots_spline[-1] - knots_spline[-neval]
    # linear spline for estimating integration
    coefs_spline <- Rmpfr::mpfrArray(NA, precBits = control$precBits,
                                     dim = c(2, neval - 1))
    coefs_spline[2, ] <- 
      (vals_spline[-1] - vals_spline[-neval]) /
      knots_diff
    coefs_spline[1, ] <- 
      vals_spline[-neval] - knots_spline[-neval] * coefs_spline[2, ]
    integral_dx <- sum(coefs_spline[1, ] * knots_spline[-1] +
                         coefs_spline[2, ] / 2 * knots_spline[-1]^2 -
                         coefs_spline[1, ] * knots_spline[-neval] - 
                         coefs_spline[2, ] / 2 * knots_spline[-neval]^2)
    # error estimation
    errors_spline <- estimate_errors(
      knots_diff, 
      c(Rmpfr::mpfr(0, precBits = control$precBits), 
        coefs_spline[2, ], 
        Rmpfr::mpfr(0, precBits = control$precBits)),
      precBits = control$precBits)
    error_dx <- sum(errors_spline)
    
    if(neval >= control$max_eval)
      break
    if(integral_dx < 0)
      stop("Negative integration values; something went wrong!")
    if(integral_dx > 0) 
      if(error_dx / abs(integral_dx) < control$rel_tol |
         error_dx < control$abs_tol)
        break
  }
  
  # modify eloga, eloga2
  integral_eloga <- sum(coefs_spline[1, ] / 2 * knots_spline[-1]^2 + 
                          coefs_spline[2, ] / 3 * knots_spline[-1]^3 -
                          coefs_spline[1, ] / 2 * knots_spline[-neval]^2 -
                          coefs_spline[2, ] / 3 * knots_spline[-neval]^3)
  error_eloga <- sum(abs(knots_diff *
                           ((vals_spline * knots_spline)[-1] - 
                              (vals_spline * knots_spline)[-neval])))
  integral_eloga2 <- sum(coefs_spline[1, ] / 3 * knots_spline[-1]^3 + 
                           coefs_spline[2, ] / 4 * knots_spline[-1]^4 -
                           coefs_spline[1, ] / 3 * knots_spline[-neval]^3 -
                           coefs_spline[2, ] / 4 * knots_spline[-neval]^4)
  error_eloga2 <- sum(abs(knots_diff *
                            ((vals_spline * knots_spline^2)[-1] - 
                               (vals_spline * knots_spline^2)[-neval])))
  
  # refit for ea
  coefs_spline <- Rmpfr::mpfrArray(NA, precBits = control$precBits,
                                   dim = c(2, neval - 1))
  coefs_spline[2, ] <- 
    ((vals_spline * exp(knots_spline))[-1] - 
       (vals_spline * exp(knots_spline))[-neval]) /
    knots_diff
  coefs_spline[1, ] <- 
    (vals_spline * exp(knots_spline))[-neval] - 
    knots_spline[-neval] * coefs_spline[2, ]
  
  integral_ea <- sum(coefs_spline[1, ] * knots_spline[-1] +
                       coefs_spline[2, ] / 2 * knots_spline[-1]^2 -
                       coefs_spline[1, ] * knots_spline[-neval] - 
                       coefs_spline[2, ] / 2 * knots_spline[-neval]^2)
  error_ea <-  sum(abs(knots_diff *
                         ((vals_spline * exp(knots_spline))[-1] - 
                            (vals_spline * exp(knots_spline))[-neval])))
  
  return(c("ea" = as.double(integral_ea / integral_dx),
           "dx" = as.double(integral_dx),
           "logLik" = as.double(log(integral_dx) - sum(log(x[x > 0]))),
           "eloga" = as.double(integral_eloga / integral_dx),
           "eloga2" = as.double(integral_eloga2 / integral_dx),
           "error_ea" = as.double(error_ea),
           "error_dx" = as.double(error_dx),
           "error_eloga" = as.double(error_eloga),
           "error_eloga2" = as.double(error_eloga2),
           "time" = as.numeric(Sys.time() - time_start, units = "secs")))
}