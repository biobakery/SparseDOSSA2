dx <- function(x, 
               pi0, mu, sigma, Omega,
               control = list()) {
  control <- do.call("integrate_control", control)
  
  ## FIXME
  log_offset <- dloga(a = x,
                      pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega)
  if(log_offset == -Inf) 
    stop("Something went wrong!")
  
  int_limits <- get_intLimits(vintegrand_dx,
                              limit_max = control$limit_max,
                              limit_min = control$limit_min,
                              step_size = control$step_size,
                              x = x, pi0 = pi0, mu = mu,
                              sigma = sigma, Omega = Omega,
                              log_offset = log_offset)

  fit_integrate <- 
    integrate(vintegrand_dx,
              lower = int_limits[1], upper = int_limits[2], 
              rel.tol = control$rel.tol, abs.tol = control$abs.tol,
              subdivisions = control$subdivisions,
              x = x, pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega,
              log_offset = log_offset)
  
  if(control$jacobian) {
    fit_integrate$value <- fit_integrate$value / prod(x[x > 0])
    fit_integrate$abs.error <- fit_integrate$abs.error / prod(x[x > 0])
  }
  
  if(control$proper) {
    fit_integrate$value <- fit_integrate$value * exp(log_offset)
    fit_integrate$abs.error <- fit_integrate$abs.error * exp(log_offset)
  }
    
  if(control$only_value)
    return(fit_integrate$value)
  else
    return(fit_integrate)
}

integrate_control <- function(limit_max = 50,
                              precBits = 100,
                              limit_min = 0.1,
                              step_size = 2,
                              rel.tol = 1e-6,
                              abs.tol = 0,
                              subdivisions = 10000,
                              jacobian = FALSE,
                              proper = TRUE,
                              only_value = TRUE) {
  list(limit_max = limit_max,
       precBits = precBits,
       limit_min = limit_min,
       step_size = step_size,
       rel.tol = rel.tol,
       abs.tol = abs.tol,
       subdivisions = subdivisions,
       jacobian = jacobian,
       proper = proper,
       only_value = only_value)
}

dloga <- function(a,
                  pi0, mu, sigma, Omega, 
                  log = TRUE) {
  u <- a_to_u(a,
              pi0 = pi0, mu = mu, sigma = sigma)
  g <- u_to_g(u = u, a = a, mu = mu, sigma = sigma) ## FIXME
  
  if(any(abs(g) == Inf)) {
    log_d <- -Inf
  } else {
    log_d <- du(g, Omega) - 
      sum((log(a[a > 0]) - mu[a > 0])^2 / (sigma[a > 0])^2 / 2) -
      log(2 * pi) / 2 * sum(a > 0) - sum(log(sigma[a > 0])) 
  }
    
  if(log)
    return(log_d)
  else
    return(exp(log_d))
}

du <- function(g, Omega, log = TRUE) {
  # Without normalizing constant (2pi)^(-2/p)!
  if(any(g == -Inf | g == Inf)) return(-Inf)
  log_d <- log_dmvnorm(S = g %*% t(g), Omega = Omega) + sum(g^2/2) + 
    log(det(Omega)) / 2
  
  if(log) 
    return(log_d)
  else
    return(exp(log_d))
}

log_dmvnorm <- function(S, Omega) {
  - sum(Omega * S) / 2
}

integrand_dx <- function(log_asum, x, 
                         pi0, mu, sigma, Omega,
                         log_offset) {
  exp(dloga(a = a(x, exp(log_asum)),
            pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega) - 
        log_offset)
}

vintegrand_dx <- Vectorize(integrand_dx, vectorize.args = "log_asum")

ea <- function(x, 
               pi0, mu, sigma, Omega,
               control) {
  control <- do.call("integrate_control", control)
  
  ## FIXME
  log_offset <- dloga(a = x,
                  pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega)
  if(log_offset == -Inf) 
    stop("Something went wrong!")
  
  int_limits <- get_intLimits(vintegrand_ea, 
                              limit_max = control$limit_max, 
                              limit_min = control$limit_min,
                              step_size = control$step_size,
                              x = x, pi0 = pi0, mu = mu, 
                              sigma = sigma, Omega = Omega,
                              log_offset = log_offset)
  
  fit_integrate <- 
    integrate(vintegrand_ea,
              subdivisions = control$subdivisions,
              rel.tol = .Machine$double.eps^0.5,
              abs.tol = 0,
              lower = int_limits[1], upper = int_limits[2],
              x = x, pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega, 
              log_offset = log_offset)  
  
  if(control$jacobian) {
    fit_integrate$value <- fit_integrate$value / prod(x[x > 0])
    fit_integrate$abs.error <- fit_integrate$abs.error / prod(x[x > 0])
  }
  
  if(control$proper) {
    fit_integrate$value <- fit_integrate$value * exp(log_offset)
    fit_integrate$abs.error <- fit_integrate$abs.error * exp(log_offset)
  }
  
  if(control$only_value)
    return(fit_integrate$value)
  else
    return(fit_integrate)
}

integrand_ea <- function(log_asum, x, 
                         pi0, mu, sigma, Omega, 
                         log_offset) {
  exp(dloga(a = a(x, exp(log_asum)),
            pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega,
            log = TRUE) - 
        log_offset + log_asum)
}

vintegrand_ea <- Vectorize(integrand_ea, 
                           vectorize.args = "log_asum")
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
# vintegrand_num_asum <- Vectorize(integrand_num_asum,
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
# vintegrand_num_asum <- Vectorize(integrand_num_asum,
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