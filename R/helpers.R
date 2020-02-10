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
      # just in case it's not pos-def 
      svd_fit$d <- abs(svd_fit$d) ## FIXME
      x_out <- svd_fit$u %*% diag(svd_fit$d) %*% t(svd_fit$u)
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

get_intLimits <- function(f, 
                          center = 0, limit_max, limit_min, step_size,
                          lower_bound = -1000, upper_bound = 1000,
                          max_try = 10,
                          ...) {
  i_try <- 0
  while(TRUE) { ## FIXME
    i_try <- i_try + 1
    if(i_try > max_try)
      stop("Could not find positive values for f!")
    
    vchange <- exp(seq(from = log(limit_max),
                       to = log(limit_min),
                       by = -log(step_size)))
    
    vlim <- c(center - vchange,
              center,
              center + rev(vchange))
    vlim <- vlim[exp(vlim) > 0] ## FIXME??
    
    vval <- f(vlim, ...)
    if(any(vval < 0))
      stop("There are negative values of f!")
    vflag <- vval > 0
    if(any(vflag))
      break
    
    step_size <- sqrt(step_size)
  }
  vindex <- which(vflag)
  if(vflag[1]) {
    warning("f is already positive at maximum lower limit!")
    lower <- lower_bound
  } else {
    lower <- vlim[min(vindex) - 1]
  }
  if(rev(vflag)[1]) {
    warning("f is already positive at maximum upper limit!")
    upper <- upper_bound
  } else {
    upper <- vlim[max(vindex) + 1]
  }
  
  return(c(lower, upper))
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