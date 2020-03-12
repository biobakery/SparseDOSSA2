logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 
  exp(x) / (1 + exp(x))

lower_tri <- function(x, warning = TRUE) {
  if(!isSymmetric(x) & warning) 
    warning("x is not symmetric!")
  
  x[lower.tri(x)]
}

upper_tri <- function(x, warning = TRUE) {
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


a_to_u <- function(a, pi0, mu, sigma) {
  to_return <-  pi0 / 2
  if(any(a > 0)) ##FIXME
    to_return[a > 0] <- 
      pi0[a > 0] + 
      pnorm(log(a[a > 0]), 
            mean = mu[a > 0], 
            sd = sigma[a > 0]) * 
      (1 - pi0[a > 0])
  
  return(to_return)
}

u_to_g <- function(u, a, mu, sigma) {
  g <- qnorm(u)
  return(g)
}

a <- function(x, asum) {
  a <- x * asum
  # In case asum is Inf
  # FIXME since this isn't supposed to happen?
  a[x == 0] <- 0  
  
  return(a)
}

get_sigmas <- function(x, eloga, eloga2, mu) {
  vapply(seq_len(ncol(x)),
         function(i_feature) {
           ind_samples <- x[, i_feature] > 0
           sqrt(mean(eloga2[ind_samples] - 
                     2 * eloga[ind_samples] * (mu[i_feature] - log(x[ind_samples, i_feature])) +
                     (mu[i_feature] - log(x[ind_samples, i_feature]))^2))
         },
         0.0)
}

# get_intLimits <- function(f, 
#                           center = 0, limit_max, limit_min, step_size,
#                           lower_bound = -1000, upper_bound = 1000,
#                           max_try = 20,
#                           ...) {
#   i_try <- 0
#   while(TRUE) { ## FIXME
#     i_try <- i_try + 1
#     if(i_try > max_try)
#       stop("Could not find positive values for f!")
#     
#     vchange <- exp(seq(from = log(limit_max),
#                        to = log(limit_min),
#                        by = -log(step_size)))
#     
#     vlim <- c(center - vchange,
#               center,
#               center + rev(vchange))
#     vlim <- vlim[exp(vlim) > 0] ## FIXME??
#     
#     vval <- f(vlim, ...)
#     if(any(vval < 0))
#       stop("There are negative values of f!")
#     vflag <- vval > 0
#     if(sum(vflag) > 1)
#       break
#     
#     step_size <- sqrt(step_size)
#   }
#   vindex <- which(vflag)
#   if(vflag[1]) {
#     warning("f is already positive at maximum lower limit!")
#     lower <- lower_bound
#   } else {
#     lower <- vlim[min(vindex) - 1]
#   }
#   if(rev(vflag)[1]) {
#     warning("f is already positive at maximum upper limit!")
#     upper <- upper_bound
#   } else {
#     upper <- vlim[max(vindex) + 1]
#   }
#   
#   return(c(lower, upper))
# }

get_intLimits <- function(x, pi0, mu, sigma, Omega,
                          n_vals = 10, step_size = 2, maxit = 10) {
  ind_nonzero <- x != 0
  logx <- log(x[ind_nonzero])
  mat_range <- cbind(mu[ind_nonzero] - logx - sqrt(2000)*sigma[ind_nonzero],
                 mu[ind_nonzero] - logx + sqrt(2000)*sigma[ind_nonzero])
  range <- c(max(mat_range[, 1]), 
             min(mat_range[, 2]))
  
  if(range[1] >= range[2])
    return(c(0, 0)) ## FIXME??

  i_iter <- 0
  while(TRUE) {
    if(i_iter + 1 > maxit)
      return(c(0, 0)) ## FIXME??
    i_iter <- i_iter + 1
    
    vlim <- seq(from = range[1], to = range[2], length.out = n_vals)
    vval <- vintegrand_dx(log_asum = vlim, x = x, pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega)
    vflag <- vval > 0
    if(vflag[1] | rev(vflag)[1])
      stop("Positive dx values at integration limits!")
    if(sum(vflag) > 1)
      return(vlim[c(min(which(vflag)) - 1, 
                    max(which(vflag)) + 1)])
    
    n_vals <- n_vals * step_size
  }
}

get_offset <- function(x, 
                       pi0, mu, sigma, Omega,
                       limit_max, limit_min, step_size) {
  
  vchange <- exp(seq(from = log(limit_max),
                     to = log(limit_min),
                     by = -log(step_size)))
  vlim <- c(-vchange, vchange)
  vlim <- vlim[exp(vlim) > 0] ## FIXME??
  
  vval <- vapply(vlim, 
                 function(i_loga) {
                   dloga(a = a(x, exp(i_loga)),
                         pi0 = pi0, mu = mu, sigma = sigma, Omega = Omega)
                 },
                 0.0)
  if(any(vval_lower < 0))
    stop("There are negative values of f!")
  vflag_lower <- vval_lower > 0
  if(any(vflag_lower)) {
    if(vflag_lower[1]) {
      warning("f is already positive at maximum lower limit!")
      lower <- lower_bound
    }
    lower <- vlim_lower[c(vflag_lower[-1], TRUE)][1]
  } else {
    warning("No positive f value for lower limits!")
    lower <- vlim_lower[length(vlim_lower)]
  }
  
  vval_upper <- f(vlim_upper, ...)
  if(any(vval_upper < 0))
    stop("There are negative values of f!")
  vflag_upper <- vval_upper > 0
  if(any(vflag_upper)) {
    if(vflag_upper[1]) {
      warning("f is already positive at maximum upper limit!")
      upper <- upper_bound
    }
    upper <- vlim_upper[c(vflag_upper[-1], TRUE)][1]
  } else {
    warning("No positive f value for upper limits!")
    upper <- vlim_upper[length(vlim_upper)]
  }
  
  return(c(lower, upper))
}

get_diff <- function(x, x_old, 
                     denom_c = 1e-5,
                     method = "abs") {
  x <- as.vector(x)
  x_old <- as.vector(x_old)
  
  abs_diff <- abs(x - x_old)
  if(method == "abs")
    return(max(abs_diff))
  if(method == "rel") {
    rel_diff <- abs_diff / (abs(x) + denom_c)
    return(max(rel_diff))
  }
}

Vectorize2 <- function(FUN, vectorize.args = arg.names, SIMPLIFY = TRUE, 
                       USE.NAMES = TRUE) 
{
  arg.names <- as.list(formals(FUN))
  arg.names[["..."]] <- NULL
  arg.names <- names(arg.names)
  vectorize.args <- as.character(vectorize.args)
  
  if(length(vectorize.args) != 1)
    stop("Can only vectorize over one argument!")
  
  if (!length(vectorize.args)) 
    return(FUN)
  if (!all(vectorize.args %in% arg.names)) 
    stop("must specify names of formal arguments for 'vectorize'")
  collisions <- arg.names %in% c("FUN", "SIMPLIFY", "USE.NAMES", 
                                 "vectorize.args")
  if (any(collisions)) 
    stop(sQuote("FUN"), " may not have argument(s) named ", 
         paste(sQuote(arg.names[collisions]), collapse = ", "))
  FUNV <- function() {
    args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    names <- if (is.null(names(args))) 
      character(length(args))
    else names(args)
    dovec <- names %in% vectorize.args
    val <- do.call("mapply", 
                   c(FUN = FUN, args[dovec], 
                     MoreArgs = list(args[!dovec]), 
                     SIMPLIFY = TRUE, USE.NAMES = USE.NAMES))
    if(!is.null(dim(args[dovec][[1]])))
      val <- array(val, dim = dim(args[dovec][[1]]))
    return(val)
  }
  formals(FUNV) <- formals(FUN)
  FUNV
}

det2 <- function(m) {
  if(!isSymmetric(m))
    stop("Matrix must be symmetric!")
  eigens <- svd(m)$d
  if(any(eigens <= 0))
    stop("Negative eigen values found for the matrix!")
  return(prod(eigens))
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
    
    ind_feature_tmp <- apply(data[ind_sample, ind_feature, drop = FALSE] > 0, 2, sum) >= k_feature
    ind_feature[ind_feature] <- ind_feature_tmp
    ind_sample_tmp <- apply(data[ind_sample, ind_feature, drop = FALSE] > 0, 1, sum) >= k_sample
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