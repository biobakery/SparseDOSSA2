dx <- function(x, 
               pi0, mu, sigma, Omega, Sigma,
               control = list(),
               log.p = FALSE) {
  control <- do.call(control_integrate, control)
  if(is.null(control$lower_loga) | is.null(control$upper_loga)) {
    limits <- get_intLimits2(x = x,
                             pi0 = pi0, mu = mu, sigma = sigma, 
                             Omega = Omega, Sigma = Sigma,
                             control = control)
    control$lower_loga <- limits[1]
    control$upper_loga <- limits[2]
  }
  
  fit_integrate <- 
    integrate2(vintegrand_dx,
               lower = control$lower_loga, upper = control$upper_loga, 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               method = control$method, 
               offset = control$offset,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  if(control$jacobian) {
    fit_integrate$integral <- fit_integrate$integral / prod(x[x > 0])
    fit_integrate$error <- fit_integrate$error / prod(x[x > 0])
  }
  
  if(log.p) 
    return(log(fit_integrate$integral))
  if(control$only_value)
    return(fit_integrate$integral)
  else
    return(fit_integrate)
}

control_integrate <- function(rel_tol = 1e-2,
                              abs_tol = 0,
                              max_eval = 50,
                              method = "cubspline",
                              step_size_limits = 2, 
                              n_vals_limits = 10,
                              maxit_limits = 4,
                              limit_tol = 1e-5,
                              offset = FALSE, 
                              lower_loga = NULL,
                              upper_loga = NULL,
                              jacobian = FALSE,
                              only_value = TRUE) {
  list(lower_loga = lower_loga,
       upper_loga = upper_loga,
       rel_tol = rel_tol,
       limit_tol = limit_tol,
       abs_tol = abs_tol,
       max_eval = max_eval,
       method = method,
       offset = offset,
       jacobian = jacobian,
       only_value = only_value,
       n_vals_limits = n_vals_limits,
       step_size_limits = step_size_limits, 
       maxit_limits = maxit_limits)
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
                       mean = rep(0, length(g)),
                       sigma = Sigma,
                       log = TRUE) +
      sum(g^2) / 2 -
      sum((log(a) - mu)^2 / (sigma)^2 / 2) - 
      sum(log(sigma)) + 
      sum(log(1 - pi0))
  } else if(!any(ind_nonzero)) {
    log_d <- 
      log(mvtnorm::pmvnorm(
        lower = -Inf, 
        upper = g, 
        mean = rep(0, length(g)),
        sigma = Sigma))
  } else {
    log_d <- 
      log(mvtnorm::pmvnorm(
        lower = -Inf, 
        upper = g[!ind_nonzero], 
        mean = (-solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE],
                       Omega[!ind_nonzero, ind_nonzero, drop = FALSE]) %*% 
                  g[ind_nonzero])[, 1],
        sigma = solve(Omega[!ind_nonzero, !ind_nonzero, drop = FALSE]))) +
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
  if(is.null(control$lower_loga) | is.null(control$upper_loga)) {
    limits <- get_intLimits2(x = x,
                             pi0 = pi0, mu = mu, sigma = sigma, 
                             Omega = Omega, Sigma = Sigma,
                             control = control)
    control$lower_loga <- limits[1]
    control$upper_loga <- limits[2]
  }
  
  fit_integrate <- 
    integrate2(vintegrand_ea,
               lower = control$lower_loga, upper = control$upper_loga, 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               method = control$method, 
               offset = control$offset,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
  
  if(control$jacobian) {
    fit_integrate$integral <- fit_integrate$integral / prod(x[x > 0])
    fit_integrate$error <- fit_integrate$error / prod(x[x > 0])
  }
  
  if(control$only_value)
    return(fit_integrate$integral)
  else
    return(fit_integrate)
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
  if(is.null(control$lower_loga) | is.null(control$upper_loga)) {
    limits <- get_intLimits2(x = x,
                             pi0 = pi0, mu = mu, sigma = sigma, 
                             Omega = Omega, Sigma = Sigma,
                             control = control)
    control$lower_loga <- limits[1]
    control$upper_loga <- limits[2]
  }
  
  fit_integrate <- 
    integrate2(vintegrand_eloga,
               lower = control$lower_loga, upper = control$upper_loga, 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               method = control$method, 
               offset = control$offset,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
 
  if(control$jacobian) {
    fit_integrate$integral <- fit_integrate$integral / prod(x[x > 0])
    fit_integrate$error <- fit_integrate$error / prod(x[x > 0])
  }
  
  if(control$only_value)
    return(fit_integrate$integral)
  else
    return(fit_integrate)
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
  if(is.null(control$lower_loga) | is.null(control$upper_loga)) {
    limits <- get_intLimits2(x = x,
                             pi0 = pi0, mu = mu, sigma = sigma, 
                             Omega = Omega, Sigma = Sigma,
                             control = control)
    control$lower_loga <- limits[1]
    control$upper_loga <- limits[2]
  }
  
  fit_integrate <- 
    integrate2(vintegrand_eloga2,
               lower = control$lower_loga, upper = control$upper_loga, 
               rel_tol = control$rel_tol, abs_tol = control$abs_tol, 
               max_eval = control$max_eval,
               method = control$method, 
               offset = control$offset,
               x = x, pi0 = pi0, mu = mu, sigma = sigma, 
               Omega = Omega, Sigma = Sigma)
 
  if(control$jacobian) {
    fit_integrate$integral <- fit_integrate$integral / prod(x[x > 0])
    fit_integrate$error <- fit_integrate$error / prod(x[x > 0])
  }
  
  if(control$only_value)
    return(fit_integrate$integral)
  else
    return(fit_integrate)
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
    n_vals = round(control$n_vals_limits / 2), 
    step_size = control$step_size_limits, 
    maxit = control$maxit_limits)
  # 
  # limits <- get_intLimits(
  #   x = data[i_sample, , drop = TRUE],
  #   pi0 = params$pi0, mu = params$mu, sigma = params$sigma, 
  #   Omega = params$Omega, Sigma = params$Sigma,
  #   control = control)
  # control <- do.call(control_integrate, control)
  
  neval <- 2
  knots_spline <- c(limits[1], limits[2])
  vals_spline <- vintegrand_dx(knots_spline, 
                               pi0 = pi0, x = x, mu = mu, sigma = sigma,
                               Omega = Omega, Sigma = Sigma)
  errors_spline <- Inf
  
  # find knots using dx
  verrors <- c()
  while(TRUE) {
    i_max_error <- which(errors_spline == max(errors_spline))[1]
    knots_spline <- c(knots_spline[seq(1, i_max_error)],
                      mean(knots_spline[c(i_max_error, i_max_error + 1)]),
                      knots_spline[seq(i_max_error + 1, neval)])
    vals_spline <- c(vals_spline[seq(1, i_max_error)],
                     vintegrand_dx(knots_spline[i_max_error + 1],
                                   pi0 = pi0, x = x, mu = mu, sigma = sigma,
                                   Omega = Omega, Sigma = Sigma),
                     vals_spline[seq(i_max_error + 1, neval)])
    
    neval <- neval + 1
    knots_diff <-  knots_spline[-1] - knots_spline[-neval]
    errors_spline <- 
      abs(knots_diff *
            (vals_spline[-1] - vals_spline[-neval]))
    
    fit_lm <- 
      lm(vals_spline ~ 
           splines::bs(knots_spline, knots = knots_spline[-c(1, neval)], degree = 1))
    coefs_spline <- 
      SplinesUtils::RegBsplineAsPiecePoly(
        fit_lm,
        "splines::bs(knots_spline, knots = knots_spline[-c(1, neval)], degree = 1)",
        shift = FALSE)$PiecePoly$coef
    # coefs_spline <-
    #   SplinesUtils::CubicInterpSplineAsPiecePoly(
    #     x = knots_spline,
    #     y = vals_spline,
    #     method = "natural")$PiecePoly$coef
    integral_dx <- sum(coefs_spline[1, ] * knots_spline[-1] +
                         coefs_spline[2, ] / 2 * knots_spline[-1]^2 -
                         coefs_spline[1, ] * knots_spline[-neval] - 
                         coefs_spline[2, ] / 2 * knots_spline[-neval]^2)
    error_dx <- sum(errors_spline)
    verrors <- c(verrors, error_dx)
    
    if(neval >= control$max_eval |
       error_dx / abs(integral_dx) < control$rel_tol |
       error_dx < control$abs_tol)
      break
  }
  
  error_dx <- sum(errors_spline)
  verrors <- c(verrors, error_dx)
  
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
  fit_lm <- 
    lm(vals_spline * exp(knots_spline) ~ 
         splines::bs(knots_spline, knots = knots_spline[-c(1, neval)], degree = 1))
  coefs_spline <- 
    SplinesUtils::RegBsplineAsPiecePoly(
      fit_lm,
      "splines::bs(knots_spline, knots = knots_spline[-c(1, neval)], degree = 1)",
      shift = FALSE)$PiecePoly$coef
  integral_ea <- sum(coefs_spline[1, ] / 3 * knots_spline[-1]^3 + 
                       coefs_spline[2, ] / 4 * knots_spline[-1]^4 -
                       coefs_spline[1, ] / 3 * knots_spline[-neval]^3 -
                       coefs_spline[2, ] / 4 * knots_spline[-neval]^4)
  error_ea <-  sum(abs(knots_diff *
                         ((vals_spline * exp(knots_spline))[-1] - 
                            (vals_spline * exp(knots_spline))[-neval])))
  
  return(c("ea" = integral_ea / integral_dx,
           "logLik" = log(integral_dx),
           "eloga" = integral_eloga / integral_dx,
           "eloga2" = integral_eloga2 / integral_dx,
           "error_ea" = error_ea,
           "error_dx" = error_dx,
           "error_eloga" = error_eloga,
           "error_eloga2" = error_eloga2,
           "time" = as.numeric(Sys.time() - time_start, units = "secs")))
}

log_intx_dc <- function(data, zero_inflation = TRUE) {
  lgamma_data <- lgamma(data + 1)
  if(!zero_inflation) {
    lgamma_sum <- lgamma(apply(data + 1, 1, sum))
    sum(lgamma_sum) - sum(lgamma_data)
  } else {
    lgamma_sum <- lgamma(apply(data, 1, function(x) sum(x[x > 0] + 1)))
    sum(lgamma_sum) - sum(lgamma_data[data > 0])
  }
}

ddirichlet <- function (x, alpha, log.p = TRUE) 
{
  
  if (!is.matrix(x)) {
    if (is.data.frame(x)) 
      x <- as.matrix(x)
  }
  
  if (!is.matrix(alpha)) 
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                    byrow = TRUE)
  if (any(dim(x) != dim(alpha))) 
    stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- ddirichlet1(x[i, ], alpha[i, ])
  if(!log.p) pd <- exp(pd)
  return(pd)
}

ddirichlet1 <- function(x, alpha) {
  logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  s <- sum((alpha - 1) * log(x))
  sum(s) - logD
}
# integrand_num_asum <- function(log_asum, x, params) {
#   asum <- exp(log_asum)
#   u <- a_to_u(a(x, asum), 
#               pi0 = params$pi0, mu = params$mu, sigma = params$sigma)
#   g <- u_to_g(u = u, a = a(x, asum), mu = params$mu, sigma = params$sigma)
#   logLik <- logLik_copula(g = g, asum = asum, x = x,
#                           mu = params$mu, sigma = params$sigma, 
#                           Omega = params$Omega)
#   
#   
#   return(exp(logLik) * asum)
# }
# 
# vintegrand_num_asum <- Vectorize2(integrand_num_asum,
#                                  vectorize.args = "log_asum")

# 
# 
# 
# 
# 
# integrand_num_asum <- function(log_asum, x, params) {
#   asum <- exp(log_asum)
#   u <- a_to_u(a(x, asum), 
#               pi0 = params$pi0, mu = params$mu, sigma = params$sigma)
#   g <- u_to_g(u = u, a = a(x, asum), mu = params$mu, sigma = params$sigma)
#   logLik <- logLik_copula(g = g, asum = asum, x = x,
#                           mu = params$mu, sigma = params$sigma, 
#                           Omega = params$Omega)
#   
#   
#   return(exp(logLik) * asum)
# }
# 
# vintegrand_num_asum <- Vectorize2(integrand_num_asum,
#                                  vectorize.args = "log_asum")
# 
# 
# pcopulasso <- function(x, 
#                        mean, sd, pi0, sigma,
#                        R = 100000) {
#   samples_a <- SparseDOSSA2:::rcopulasso(n = R,
#                                          mean = params$mu,
#                                          sd = params$sigma,
#                                          pi0 = params$pi0,
#                                          sigma = solve(params$Omega))
#   samples_asum <- vapply(seq_len(R), function(r) sum(samples_a[r, ]), 
#                          0.0)
#   
#   for(i in seq_len(nrow(x))) {
#     i_rx <- rx(c = x[i, ], C = sum(x[i, ]), R = R)
#     samples_x <- i_rx$samples
#     samples_xsum <- vapply(seq_len(R), function(r) sum(samples_a[r, ]), 
#                            0.0)
#   }
#   l_rx <- lapply(seq_len(nrow(x)), 
#                  function(i) {
#                    rx(c = x[i, ], C = sum(x[i, ]), R = R)
#                  })
#   
#   sum(log(vapply(seq_len(nrow(x)), 
#                  function(i)
#                    mean(vapply(seq_along(R), function(r) {
#                      dmultinom(x = x[i, ], prob = samples_X[r, ])
#                    }, 0.0)),
#                  0.0)))
# }
# 
# rx <- function(c, C, R) {
#   x <- (c + 0.5) / (C + 0.5)
#   mu <- log(x)
#   sigma <- sqrt((1 - x) / (c + 0.5))
#   
#   samples <- exp(rnorm(n = R * length(c), mean = mu, sd = sigma)) * 
#     rbinom(n = R * length(c), size = 1, prob = 0.5)
#   samples <- t(matrix(samples, nrow = length(c)))
#   params <- list(mu = mu, sigma = sigma)
#   
#   return(list(samples = samples, params = params))
# }