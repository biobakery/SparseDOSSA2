dx <- function(x, 
               pi0, mu, sigma, Omega, Sigma,
               control = list(),
               log.p = FALSE) {
  control <- do.call(control_integrate, control)
  
  int_limits <- get_intLimits(x = x, pi0 = pi0, mu = mu, sigma = sigma, 
                              Omega = Omega, Sigma = Sigma,
                              maxit = control$maxit_getLimits)
  
  fit_integrate <- 
    cubature::cubintegrate(vintegrand_dx,
                           lower = int_limits[1], upper = int_limits[2], 
                           relTol = control$rel_tol, absTol = control$abs_tol,
                           method = control$method, maxEval = control$max_eval,
                           nVec = 2,
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

control_integrate <- function(limit_max = 50,
                              limit_min = 1e-5,
                              step_size = 2,
                              maxit_getLimits = 10,
                              rel_tol = 1e-05,
                              abs_tol = 0,
                              max_eval = 1e6,
                              method = "hcubature",
                              jacobian = FALSE,
                              only_value = TRUE) {
  list(limit_max = limit_max,
       limit_min = limit_min,
       step_size = step_size,
       maxit_getLimits = maxit_getLimits,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       max_eval = max_eval,
       method = method,
       jacobian = jacobian,
       only_value = only_value)
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
      sum((log(a[ind_nonzero]) - mu[ind_nonzero])^2 / (sigma[ind_nonzero])^2 / 2) - 
      sum(log(sigma[ind_nonzero])) + 
      sum(log(1 - pi0[ind_nonzero]))
  }
  
  if(log.p)
    return(log_d)
  else
    return(exp(log_d))
}

# dloga_old <- function(a,
#                       pi0, mu, sigma, Omega, 
#                       log.p = TRUE) {
#   u <- a_to_u_old(a,
#                   pi0 = pi0, mu = mu, sigma = sigma)
#   g <- qnorm(u)
#   
#   if(any(abs(g) == Inf)) {
#     log_d <- -Inf
#   } else {
#     log_d <- du(g, Omega) - 
#       sum((log(a[a > 0]) - mu[a > 0])^2 / (sigma[a > 0])^2 / 2) -
#       log(2 * pi) / 2 * sum(a > 0) - sum(log(sigma[a > 0])) 
#   }
#   
#   if(log.p)
#     return(log_d)
#   else
#     return(exp(log_d))
# }
# 
# du <- function(g, Omega, log.p = TRUE) {
#   # Without normalizing constant (2pi)^(-2/p)!
#   if(any(g == -Inf | g == Inf)) return(-Inf)
#   log_d <- log_dmvnorm(S = g %*% t(g), Omega = Omega) + sum(g^2)/2 + 
#     log(det(Omega)) / 2
#   
#   if(log.p) 
#     return(log_d)
#   else
#     return(exp(log_d))
# }
# 
# log_dmvnorm <- function(S, Omega) {
#   - sum(Omega * S) / 2
# }

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
               control) {
  control <- do.call(control_integrate, control)
  
  int_limits <- get_intLimits(x = x, pi0 = pi0, mu = mu, sigma = sigma, 
                              Omega = Omega, Sigma = Sigma,
                              maxit = control$maxit_getLimits)
  
  fit_integrate <- 
    cubature::cubintegrate(vintegrand_ea,
                           lower = int_limits[1], upper = int_limits[2], 
                           relTol = control$rel_tol, absTol = control$abs_tol,
                           method = control$method, maxEval = control$max_eval,
                           nVec = 2,
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
                  control) {
  control <- do.call(control_integrate, control)
  
  int_limits <- get_intLimits(x = x, pi0 = pi0, mu = mu, sigma = sigma, 
                              Omega = Omega, Sigma = Sigma,
                              maxit = control$maxit_getLimits)
  
  fit_integrate <- 
    cubature::cubintegrate(vintegrand_eloga,
                           lower = int_limits[1], upper = int_limits[2], 
                           relTol = control$rel_tol, absTol = control$abs_tol,
                           method = control$method, maxEval = control$max_eval,
                           nVec = 2,
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
                   control) {
  control <- do.call(control_integrate, control)
  
  int_limits <- get_intLimits(x = x, pi0 = pi0, mu = mu, sigma = sigma, 
                              Omega = Omega, Sigma = Sigma,
                              maxit = control$maxit_getLimits)
  
  fit_integrate <- 
    cubature::cubintegrate(vintegrand_eloga2,
                           lower = int_limits[1], upper = int_limits[2], 
                           relTol = control$rel_tol, absTol = control$abs_tol,
                           method = control$method, maxEval = control$max_eval,
                           nVec = 2,
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