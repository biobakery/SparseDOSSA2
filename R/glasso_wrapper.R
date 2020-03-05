glasso_wrapper <- function(S, lambda, source = "huge",
                           simplify = FALSE,
                           symm = TRUE,
                           corr = TRUE,
                           threshold_zero = 1e-16,
                           debug_file = NULL) {
  Omega <- diag(1/(diag(S) + lambda))
  
  if(simplify) 
    z <- which(rowSums(abs(S) > lambda) > 1)
  else 
    z <- seq_len(nrow(S))
  q <- length(z)
  if (q > 0) {
    if(source == "huge") {
      out.glasso <- huge::huge.glasso(x = S[z, z, drop = FALSE],
                                      lambda = lambda,
                                      verbose = FALSE)
      Omega[z, z] <- out.glasso$icov[[1]]
    }
    if(source == "hugec") {
      out.glasso <- .Call("_huge_hugeglasso",
                          S[z, z, drop = FALSE],
                          lambda,
                          FALSE,
                          FALSE,
                          FALSE,
                          PACKAGE = "huge")
      Omega[z, z] <- out.glasso$icov[[1]]
    }
    if(source == "glasso") {
      out.glasso <- glasso::glasso(s = S[z, z, drop = FALSE],
                                   rho = lambda)
      Omega[z, z] <- out.glasso$wi
    }
  }
  
  if(!is.null(debug_file))
    save(Omega, file = debug_file)

  if(any(is.na(Omega))) {
    print("Missing values in Omega estimation! (lambda to small?)") # FIXME
    stop("Missing values in Omega estimation! (lambda to small?)")
  }
  
  if(symm) {
    Omega <- enforce_symm(Omega, method = "svd")
  }
  if(corr) {
    Omega <- enforce_corr(Omega)
  }
  Omega <- threshold_matrix(Omega, threshold_zero = threshold_zero)
  
  return(Omega)
}
