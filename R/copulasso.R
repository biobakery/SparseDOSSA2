copulasso <- function(data, marginals,
                      lambda, 
                      control = list()) {
  control <- do.call(control_copulasso, control)
  if(ncol(data) != nrow(marginals))
    stop("Dimension of data and marginals do not agree!")
  
  if(is.null(lambda)) {
    Omega <- 
      Sigma <- 
      Corr_star <- diag(1, ncol(data))
  } else {
    data_g <- vapply(seq_len(ncol(data)),
                     function(i_feature)
                       a_to_g(a = data[, i_feature],
                              pi0 = marginals[i_feature, "pi0"],
                              mu = marginals[i_feature, "mu"],
                              sigma = marginals[i_feature, "sigma"]),
                     rep(0.0, nrow(data)))
    Corr_star <- cor(data_g, method = "pearson")
    S <- get_s(cor_data = Corr_star,
               pi0 = marginals[, "pi0"],
               glim = marginals[, "glim"],
               g0 = marginals[, "g0"],
               sigmaMod = marginals[, "sigmaMod"])
    Omega <- diag(1/(diag(S) + lambda))
    
    if(control$simplify) 
      z <- which(rowSums(abs(S) > lambda) > 1)
    else 
      z <- seq_len(nrow(S))
    q <- length(z)
    if (q > 0) {
      if(control$penalize_method == "huge_conditioned") {
        S_conditioned <- condition_ridge(S[z, z, drop = FALSE],
                                         lambda = 1e-6,
                                         method = "ridge1")
        out.glasso <- huge::huge.glasso(x = S_conditioned,
                                        lambda = lambda,
                                        verbose = FALSE)
        Omega[z, z] <- out.glasso$icov[[1]]
      }
      if(control$penalize_method %in% c("ridge1", "ridge2")) {
        out.glasso <- solve_ridge(S[z, z, drop = FALSE],
                                  lambda,
                                  method = penalize_method)
        Omega[z, z] <- out.glasso
      }
    }
  }
  
  if(!is.null(control$debug_dir))
    save(Omega, file = paste0(control$debug_dir, "/debug_copulasso.RData"))

  if(any(is.na(Omega))) {
    # warning("Missing values in Omega estimation! (lambda to small?)") # FIXME
    Omega <- diag(rep(1, nrow(Omega)))
    return(list(Omega = Omega,
                copulasso_code = 1))
  }
  
  if(control$symm) {
    Omega <- enforce_symm(Omega, method = "svd")
  }
  if(control$corr) {
    Omega <- enforce_corr(Omega)
  }
  
  Omega <- threshold_matrix(Omega, 
                            threshold_zero = control$threshold_zero)
  Sigma <- threshold_matrix(solve(Omega), 
                            threshold_zero = control$threshold_zero)
  
  return(list(Omega = Omega,
              Sigma = Sigma,
              Corr_star = Corr_star,
              copulasso_code = 0))
}

control_copulasso <- function(penalize_method = "huge_conditioned",
                              threshold_zero = 1e-16,
                              simplify = FALSE,
                              symm = TRUE,
                              corr = TRUE,
                              debug_dir = NULL) {
  list(penalize_method = penalize_method,
       threshold_zero = threshold_zero,
       simplify = simplify,
       symm = symm,
       corr = corr, 
       debug_dir = debug_dir)
}

Rho <- function(rho_p, 
                pi0_1, pi0_2, 
                glim_1, glim_2,
                g0_1, g0_2,
                sigmaMod_1, sigmaMod_2) {
  mat_cor <- matrix(c(1, rho_p, rho_p, 1),
                    2, 2)
  part1 <- g0_1 * g0_2 * 
    mvtnorm::pmvnorm(upper = c(glim_1, glim_2),
                     corr = mat_cor)
  part2 <- g0_1 *
    tmvtnorm::mtmvnorm(sigma = mat_cor,
                       lower = c(-Inf, glim_2),
                       upper = c(glim_1, Inf),
                       doComputeVariance = FALSE)$tmean[2] *
    mvtnorm::pmvnorm(lower = c(-Inf, glim_2),
                     upper = c(glim_1, Inf),
                     corr = mat_cor)
  part3 <- 
    g0_2 *
    tmvtnorm::mtmvnorm(sigma = mat_cor,
                       lower = c(glim_1, -Inf),
                       upper = c(Inf, glim_2),
                       doComputeVariance = FALSE)$tmean[1] *
    mvtnorm::pmvnorm(lower = c(glim_1, -Inf),
                     upper = c(Inf, glim_2),
                     corr = mat_cor)
  part4 <- {
    fit_mtmvnorm <- tmvtnorm::mtmvnorm(sigma = mat_cor,
                                       lower = c(glim_1, glim_2),
                                       doComputeVariance = TRUE)
    if (any(is.na(fit_mtmvnorm$tmean)))
      0
    else
      (prod(fit_mtmvnorm$tmean) + fit_mtmvnorm$tvar[1, 2]) *
      mvtnorm::pmvnorm(lower = c(glim_1, glim_2),
                       corr = mat_cor)
  }
  
  return((part1 + part2 + part3 + part4) / sigmaMod_1 / sigmaMod_2)
}

vRho <- Vectorize(Rho, vectorize.args = "rho_p")

iRho <- function(rho_s, 
                 pi0_1, pi0_2, 
                 glim_1, glim_2,
                 g0_1, g0_2,
                 sigmaMod_1, sigmaMod_2) {
  f_lim <- vRho(c(-0.99, 0.99),
                pi0_1 = pi0_1, pi0_2 = pi0_2, 
                glim_1 = glim_1, glim_2 = glim_2,
                g0_1 = g0_1, g0_2 = g0_2,
                sigmaMod_1 = sigmaMod_1, sigmaMod_2 = sigmaMod_2)
  if(rho_s <= f_lim[1])
    return(-0.99) ## FIXME
  if(rho_s >= f_lim[2])
    return(0.99)
  
  uniroot(f = function(x) 
    vRho(x, 
         pi0_1 = pi0_1, pi0_2 = pi0_2, 
         glim_1 = glim_1, glim_2 = glim_2,
         g0_1 = g0_1, g0_2 = g0_2,
         sigmaMod_1 = sigmaMod_1, sigmaMod_2 = sigmaMod_2) - 
      rho_s, 
    interval = c(-0.99, 0.99),
    f.lower = f_lim[1],
    f.upper = f_lim[2],
    maxiter = 100)$root
}

get_s <- function(cor_data, pi0, glim, g0, sigmaMod) {
  df_index <- expand.grid(feature2 = seq_len(nrow(cor_data)),
                          feature1 = seq_len(nrow(cor_data)))
  df_index <- subset(df_index, feature1 < feature2)
  df_index$i_combo <- seq_len(nrow(df_index))
  
  ss <- 
    future.apply::future_vapply(
      df_index$i_combo,
      function(ii_combo) {
        ind_feature1 <- df_index[ii_combo, ]$feature1
        ind_feature2 <- df_index[ii_combo, ]$feature2
        
        iRho(rho_s = cor_data[ind_feature1, ind_feature2],
             pi0_1 = pi0[ind_feature1], pi0_2 = pi0[ind_feature2],
             glim_1 = glim[ind_feature1], glim_2 = glim[ind_feature2],
             g0_1 = g0[ind_feature1], g0_2 = g0[ind_feature2],
             sigmaMod_1 = sigmaMod[ind_feature1], 
             sigmaMod_2 = sigmaMod[ind_feature2])
      },
      0.0)
  
  to_return <- cor_data
  to_return[lower.tri(to_return)] <- ss
  to_return[upper.tri(to_return)] <- upper_tri(t(to_return), warning = FALSE)
  
  return(to_return)
}

solve_ridge <- function(S, lambda, method = "ridge1") {
  svd_fit <- svd(S)
  if(method == "ridge1")
    eigen_inv <- (sqrt(svd_fit$d^2 + 4 * lambda) - svd_fit$d) / 2 / lambda
  if(method == "ridge2")
    eigen_inv <- 1 / (svd_fit$d + lambda)
  
  return(svd_fit$u %*% diag(eigen_inv) %*% t(svd_fit$u))
}

condition_ridge <- function(S, lambda, method = "ridge1") {
  svd_fit <- svd(S)
  if(method == "ridge1")
    eigen_cond <- 2 * lambda / (sqrt(svd_fit$d^2 + 4 * lambda) - svd_fit$d)
  if(method == "ridge2")
    eigen_cond <- 1 / (svd_fit$d + lambda)
  
  return(svd_fit$u %*% diag(eigen_cond) %*% t(svd_fit$u))
}