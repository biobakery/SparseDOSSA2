logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 
  exp(x) / (1 + exp(x))

TSS <- function(x, correct = FALSE) {
  if(any(x < 0))
    stop("Negative x values are not accepted!")
  if(all(x == 0)) return(x)
  # this is a special case where only one feature is present
  if(sum(x > 0) == 1 & correct) {
    x[x > 0] <- 1 - 0.5 / length(x)
    return(x)
  }
  return(x / sum(x))
}

lower_tri <- function(x, warning = FALSE) {
  if(!isSymmetric(x) & warning) 
    warning("x is not symmetric!")
  
  x[lower.tri(x)]
}

upper_tri <- function(x, warning = FALSE) {
  if(!isSymmetric(x) & warning) 
    warning("x is not symmetric!")
  
  x[upper.tri(x)]
}

enforce_symm <- function(x, method = "svd") {
  if(nrow(x) != ncol(x)) 
    stop("x does not appear to be a covariance matrix!")
  x_out <- x
  if(!isSymmetric(x_out)) {
    if(method == "average") {
      lower_averaged <- (lower_tri(x_out, warning = FALSE) + 
                           lower_tri(t(x_out), warning = FALSE)) / 2
      x_out[lower.tri(x_out)] <- lower_averaged
      x_out[upper.tri(x_out)] <- 
        t(x_out)[upper.tri(x_out)]
    }
    if(method == "lower")
      x_out[upper.tri(x_out)] <- upper_tri(t(x_out), warning = FALSE)
    if(method == "upper")
      x_out[lower.tri(x_out)] <- lower_tri(t(x_out), warning = FALSE)
    if(method == "svd") {
      svd_fit <- svd(x)
      if(Matrix::rankMatrix(svd_fit$u) < nrow(x)) 
        # Sometimes svd run into weird cases where the matrix is not singular but the 
        # eigen vector matrix is
        svd_fit <- svd(x + diag(rep(1e-16, nrow(x))))
      
      # In case it's not pos-def 
      svd_fit$d <- abs(svd_fit$d) ## FIXME
      
      x_out <- svd_fit$u %*% diag(svd_fit$d) %*% t(svd_fit$u)
      
      # I don't know how this could happen, but even after this x_out can still be not symmetric!
      ## FIXME
      x_out[upper.tri(x_out)] <- upper_tri(t(x_out), warning = FALSE)
    }
  }
  
  return(x_out)
}

enforce_corr <- function(x) {
  if(nrow(x) != ncol(x)) 
    stop("x does not appear to be a covariance matrix!")

  diag_sigma <- diag(sqrt(abs(diag(solve(x)))))
  x_out <- diag_sigma %*% x %*% diag_sigma
  
  return(x_out)
}

a_to_u <- function(a, pi0, mu, sigma,
                   half_pi0 = FALSE) {
  if(any(pi0 == 1) | any(pi0 == 0))
    stop("zero or one valued pi0 is not supported!")
  
  if(half_pi0)
    to_return <-  pi0 / 2
  else
    to_return <-  pi0

  ind_nonzero <- a > 0
  if(any(ind_nonzero)) ##FIXME
    to_return[ind_nonzero] <-
      pi0[ind_nonzero] +
      pnorm(log(a[ind_nonzero]),
            mean = mu[ind_nonzero],
            sd = sigma[ind_nonzero]) *
      (1 - pi0[ind_nonzero])

  return(to_return)
}

a_to_u_marg <- function(a, pi0, mu, sigma,
                    half_pi0 = FALSE,
                    u_tol =  1e-5) {
  if(half_pi0)
    to_return <-  rep(pi0 / 2,
                      length = length(a))
  else
    to_return <-  rep(pi0,
                      length = length(a))
  
  ind_nonzero <- a > 0
  if(any(ind_nonzero)) ##FIXME
    to_return[ind_nonzero] <-
    pi0 +
    pnorm(log(a[ind_nonzero]),
          mean = mu,
          sd = sigma) *
    (1 - u_tol) *
    (1 - pi0)
  
  return(to_return)
}

a_to_g <- function(a, pi0, mu, sigma) {
  if(pi0 == 1 | pi0 == 0)
    stop("zero or one valued pi0 is not supported!")
  ind_nonzero <- a > 0
  
  to_return <- u <- a_to_u_marg(a, pi0, mu, sigma)
  to_return[ind_nonzero] <- qnorm(u[ind_nonzero])
  to_return[!ind_nonzero] <- 
    -mean(to_return[ind_nonzero]) * 
    sum(ind_nonzero) / sum(1 - ind_nonzero)
  
  return(to_return)
}

a <- function(x, asum) {
  a <- x * asum
  a[x == 0] <- 0  
  
  return(a)
}

make_CVfolds <- function(n, K) {
  cut(sample.int(n), breaks = K, labels = FALSE)
}

filter_data <- function(data, 
                        k_feature = 2, k_sample = 1,
                        maxit = 3) {
  i_iter <- 0
  ind_feature <- rep(TRUE, ncol(data))
  ind_sample <- rep(TRUE, nrow(data))
  
  while(TRUE) {
    if (i_iter + 1 > maxit) 
      stop("Max iteration reached!")
    i_iter <- i_iter + 1
    
    ind_feature_tmp <- 
      apply(data[ind_sample, ind_feature, drop = FALSE] > 0, 2, sum) >= 
      k_feature
    ind_feature[ind_feature] <- ind_feature_tmp
    ind_sample_tmp <- 
      apply(data[ind_sample, ind_feature, drop = FALSE] > 0, 1, sum) >= 
      k_sample
    ind_sample[ind_sample] <- ind_sample_tmp
    
    if (all(ind_feature_tmp) & all(ind_sample_tmp)) 
      return(list(ind_feature = ind_feature,
                  ind_sample = ind_sample))
  }
}

fill_estimates_CV <- function(params_CV, params_full, ind_feature) {
  params_return <- params_CV
  
  params_return[c("pi0", "mu", "sigma")] <- params_full[c("pi0", "mu", "sigma")]
  for(param in c("pi0", "mu", "sigma"))
    params_return[[param]][ind_feature] <- params_CV[[param]]
  
  params_return$Sigma <- 
    params_return$Omega <- 
    diag(rep(1, length(ind_feature)))
  params_return[["Sigma"]][ind_feature, ind_feature] <- params_CV[["Sigma"]]
  params_return[["Omega"]][ind_feature, ind_feature] <- threshold_matrix(solve(params_CV[["Sigma"]]))
  
  return(params_return)
}

threshold_matrix <- function(x, threshold_zero = 1e-16) {
  ## FIXME
  if(!is.null(threshold_zero)) {
    diag_x <- diag(x)
    x[abs(x) < threshold_zero] <- 0
    diag(x) <- diag_x
  }
  return(x)
}