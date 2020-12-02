filter_params <- function(params, max_pi0 = 0.95) {
  ind <- params$pi0 < max_pi0
  return(list(pi0 = params$pi0[ind],
              mu = params$mu[ind],
              sigma = params$sigma[ind],
              Omega = threshold_matrix(solve(params$Sigma[ind, ind])),
              Sigma = params$Sigma[ind, ind]))
}

params_a_to_x <- function(params_a, n = 100000) {
  samples_a <- rcopulasso(n = round(n * 1.5), # FIXME
                          pi0 = params_a$pi0,
                          mu = params_a$mu,
                          sigma = params_a$sigma,
                          Omega = params_a$Omega)
  samples_a <- samples_a[!apply(samples_a == 0, 1, all), ][1:n, ]
  samples_x <- t(apply(samples_a, 1, function(x) x / sum(x)))
  mat_marginals <- get_marginals(samples_x)
  return(list(pi0 = mat_marginals[, "pi0"],
              mu = mat_marginals[, "mu"],
              sigma = mat_marginals[, "sigma"],
              Corr = cor(samples_x, method = "spear")))
}