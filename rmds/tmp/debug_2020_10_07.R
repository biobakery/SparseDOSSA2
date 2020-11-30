rm(list = ls())
library(magrittr)
library(ggplot2)
smar::sourceDir("R/")

load("../../sparsedossa_paper/data/physeqs/vaginal.RData")
x_samples <- physeq_vaginal_bl %>% 
  smar::otu_table2()
lambda <- 1^(-0.5)

load("../../sparsedossa_paper/results/fitted/SparseDOSSA2/4/K_0/lambda_1/debug_EM.RData")
l_debug$ll_params %>% 
  purrr::map_dbl("logLik") %>% 
  plot()
l_debug$ll_params %>% 
  purrr::map_dbl(~.x$mu[2]) %>% 
  plot()


# is the gradient descent correct?
data <- t(x_samples)
control <- do.call(control_fit, list())

e_asums <- l_debug$ll_easums[[10]]
params <- l_debug$ll_params[[9]]
params_new <- l_debug$ll_params[[10]]

future::plan(future::multisession())
e_asums_new <- future.apply::future_vapply(
  seq_len(nrow(data)),
  function(i_sample)
    get_es(x = data[i_sample, , drop = TRUE],
           pi0 = params$pi0, mu = params_new$mu, sigma = params_new$sigma, 
           Omega = params$Omega, Sigma = params$Sigma,
           control = control$control_numint),
  rep(0.0, 10)
) %>% t()

sum(e_asums[, "logLik"])
sum(e_asums_new[, "logLik"])

# so something went wrong during the gradient descent. confirm that.
logLik_marginal <- function(mutilde, sigma2tilde,
                            mu, sigma2) {
  sum(- 1 / 2 / sigma2 * (sigma2tilde - 2 * mutilde * mu + mu^2) -
        log(sigma2) / 2)
}

logLik_marginal2 <- function(mutilde, sigma2tilde,
                             n_nonzero,
                            sigma2) {
  sum(-log(sigma2) * n_nonzero  - (sigma2tilde - mutilde^2) / sigma2 * n_nonzero) / 2 -
    1 / 2 * sum(mutilde)^2 / sum(sigma2 / n_nonzero)
}

logLik_marginal2.5 <- function(mutilde, sigma2tilde,
                             sigma2) {
  sum(-log(sigma2)  - (sigma2tilde - mutilde^2) / sigma2) / 2 -
    1 / 2 * sum(mutilde)^2 / sum(sigma2)
}

logLik_marginal3 <- function(mutilde, sigma2tilde,
                            mu, sigma2) {
  sum(- 1 / 2 / sigma2 * ((mu - mutilde)^2 + sigma2tilde - mutilde^2) -
        log(sigma2) / 2)
}

sum(- 1 / 2 * (params_new$mu - mutilde)^2 / params_new$sigma^2)
- 1 / 2 * sum(mutilde)^2 / sum(params_new$sigma^2)

solve_mu


plot(muhat - solve_mu(mutilde, sigma2hat))

mutilde <- get_mutilde(data = data,
                       e_asums = e_asums)
sigma2tilde <- get_sigma2tilde(data = data,
                               e_asums = e_asums)
n_nonzero <- vapply(seq_len(ncol(data)),
                    function(j) sum(data[, j] > 0),
                    0)
sigma2hat <- solve_sigma2(mutilde = mutilde, 
                          sigma2tilde = sigma2tilde, 
                          n_nonzero = n_nonzero,
                          control = control, 
                          maxiter = 1000)
muhat <- solve_mu(mutilde = mutilde, 
                  sigma2hat = sigma2hat)

test <- function(x, y) {
  sapply(seq_along(x),
         function(i) logLik_marginal2(c(1, -1), sigma2tilde[1:2], n_nonzero[1:2], c(x[i], y[i])))
  
}

x <- seq(100, 10000, length.out = 100)
y <- seq(100, 10000, length.out = 100)
z <- outer(x, y, test)

persp(x, y, z)

logLik_marginal(mutilde, sigma2tilde, params_new$mu, params_new$sigma^2)
logLik_marginal(mutilde, sigma2tilde, params$mu, params$sigma^2)
logLik_marginal2(mutilde, sigma2tilde, params_new$sigma^2)
logLik_marginal2(mutilde, sigma2tilde, params$sigma^2)
logLik_marginal3(mutilde, sigma2tilde, params_new$mu, params_new$sigma^2)
logLik_marginal3(mutilde, sigma2tilde, params$mu, params$sigma^2)
logLik_marginal(mutilde, sigma2tilde, muhat, sigma2hat)
logLik_marginal2(mutilde, sigma2tilde, sigma2hat)


# initialization, filtering
if(is.null(lambda) & any(lambda <= 0))
  stop("lambda must be a positive value!")
control <- do.call(control_fit, control)
debug_copulasso_file <- NULL
if(!is.null(control$debug_dir)) {
  dir.create(control$debug_dir, recursive = TRUE)
  control$debug_dir <- normalizePath(control$debug_dir)
}

l_filtering <- filter_data(data)
data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]

# initialize EM using relative abundances
feature_param <- fit_featureParam(data)
feature_param[, "mu"] <- feature_param[, "mu"] - mean(feature_param[, "mu"])
fit_copulasso <- copulasso(
  data = data, 
  marginals = feature_param,
  lambda = lambda,
  control = c(control$control_copulasso,
              list(debug_dir = control$debug_dir))
)
params <- list(pi0 = feature_param[ ,"pi0"],
               mu = feature_param[, "mu"],
               sigma = feature_param[, "sigma"],
               Sigma = fit_copulasso$Sigma,
               Omega = fit_copulasso$Omega,
               Corr_star = fit_copulasso$Corr_star,
               diff = rep(NA_real_, 4),
               logLik = -Inf,
               time = Sys.time())
if(fit_copulasso$copulasso_code != 0) {
  warning("Missing values in Omega estimation! (lambda to small?)")
  return(list(lambda = lambda,
              fit = params,
              convergence = list(converge = converge,
                                 converge_code = 4,
                                 n_iter = i_iter)))
}



## M step
sigma2hat <- solve_sigma2(mutilde = mutilde, 
                          sigma2tilde = sigma2tilde, 
                          control = control, 
                          maxiter = 1000)
  muhat <- solve_mu(mutilde = mutilde, 
                    sigma2hat = sigma2hat)
  feature_param[, "sigma"] <- sqrt(sigma2hat)
  feature_param[, "mu"] <- muhat
  fit_copulasso <- copulasso(data = data * exp(e_asums[, "eloga"]), 
                             marginals = feature_param,
                             lambda = lambda,
                             control = control$control_copulasso)
  if(fit_copulasso$copulasso_code != 0) {
    warning("Missing values in Omega estimation! (lambda to small?)")
    converge_code <- 4
    break
  }
  params_new <- list(pi0 = feature_param[ ,"pi0"],
                     mu = feature_param[ ,"mu"],
                     sigma = feature_param[ ,"sigma"],
                     Sigma = fit_copulasso$Sigma,
                     Omega = fit_copulasso$Omega,
                     Corr_star = fit_copulasso$Corr_star,
                     logLik = mean(e_asums[, "logLik"]))
  diff_abs <- get_diff(params_new[["logLik"]], params[["logLik"]], 
                       denom_c = control$abs_tol, method = "abs")
  diff_rel <- get_diff(params_new[["logLik"]], params[["logLik"]], 
                       denom_c = control$abs_tol, method = "rel")
  params <- c(params_new,
              list(diff = c(diff_abs, diff_rel),
                   time = Sys.time()))
  ll_params[[i_iter]] <- params
  
  if(!is.null(control$debug_dir)) {
    l_debug <- list(ll_easums = ll_easums, 
                    ll_params = ll_params, 
                    l_filtering = l_filtering)
    save(l_debug,
         file = paste0(control$debug_dir,
                       "/debug_EM.RData"))
  }
  
  if(max(diff_abs) < control$abs_tol & max(diff_rel) < control$rel_tol) {
    converge <- TRUE
    converge_code <- 0
    break
  }
  if(i_iter + 1 > control$maxit) {
    warning("Maximum EM iteration reached!")
    converge_code <- 1
    break
  }
}