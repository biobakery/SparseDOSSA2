copulasso <- function(data, lambda, 
                      penalize_method = "huge",
                      simplify = FALSE,
                      symm = TRUE,
                      corr = TRUE,
                      threshold_zero = 1e-16,
                      debug_file = NULL) {
  S <- get_s(data = data)
  Omega <- diag(1/(diag(S) + lambda))
  
  if(simplify) 
    z <- which(rowSums(abs(S) > lambda) > 1)
  else 
    z <- seq_len(nrow(S))
  q <- length(z)
  if (q > 0) {
    if(penalize_method == "huge") {
      out.glasso <- huge::huge.glasso(x = S[z, z, drop = FALSE],
                                      lambda = lambda,
                                      verbose = FALSE)
      Omega[z, z] <- out.glasso$icov[[1]]
    }
    if(penalize_method == "hugec") {
      out.glasso <- .Call("_huge_hugeglasso",
                          S[z, z, drop = FALSE],
                          lambda,
                          FALSE,
                          FALSE,
                          FALSE,
                          PACKAGE = "huge")
      Omega[z, z] <- out.glasso$icov[[1]]
    }
    if(penalize_method == "glasso") {
      out.glasso <- glasso::glasso(s = S[z, z, drop = FALSE],
                                   rho = lambda)
      Omega[z, z] <- out.glasso$wi
    }
    if(penalize_method %in% c("ridge1", "ridge2")) {
      out.glasso <- solve_ridge(S[z, z, drop = FALSE],
                                lambda,
                                method = penalize_method)
      Omega[z, z] <- out.glasso
    }
  }
  
  if(!is.null(debug_file))
    save(Omega, file = debug_file)

  if(any(is.na(Omega))) {
    # warning("Missing values in Omega estimation! (lambda to small?)") # FIXME
    Omega <- diag(rep(1, nrow(Omega)))
    return(list(Omega = Omega,
                copulasso_code = 1))
  }
  
  if(symm) {
    Omega <- enforce_symm(Omega, method = "svd")
  }
  if(corr) {
    Omega <- enforce_corr(Omega)
  }
  Omega <- threshold_matrix(Omega, threshold_zero = threshold_zero)
  
  return(list(Omega = Omega,
              copulasso_code = 0))
}

iRho <- function(rho_s) sinpi(rho_s/6) * 2

Rho <- function(rho_p) asin(rho_p / 2) / pi * 6

get_s <- function(data, method = "spearman",
                  random = TRUE, sim = FALSE,
                  R = 1000) {
  s_s <- cor(x = data, method = method)
  iRho(s_s)
}

negLogLik_mvn <- function(S, Omega) { ##FIXME
  -log(det(Omega)) + sum(Omega * S)
}

solve_ridge <- function(S, lambda, method = "ridge1") {
  svd_fit <- svd(S)
  if(method == "ridge1")
    eigen_inv <- (sqrt(svd_fit$d^2 + 4 * lambda) - svd_fit$d) / 2 / lambda
  if(method == "ridge2")
    eigen_inv <- 1 / (svd_fit$d + lambda)
  
  return(svd_fit$u %*% diag(eigen_inv) %*% t(svd_fit$u))
}