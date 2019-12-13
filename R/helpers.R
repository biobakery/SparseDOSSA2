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

enforce_symm <- function(x, method = "upper") {
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
  }
  
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
  a <- x
  a[x > 0] <- x[x > 0] * asum # in case asum is Inf
  
  return(a)
}

get_sigmas <- function(data, mu) {
  vapply(seq_len(ncol(data)),
         function(i_sample) {
           sqrt(mean((log((data[, i_sample])[data[, i_sample] > 0]) - 
                        mu[i_sample])^2))
         },
         0.0)
}

get_intLimits <- function(f, 
                          limit_max, limit_min, step_size,
                          ...) {
  vlim_lower <- -exp(seq(from = log(limit_max),
                         to = log(limit_min),
                         by = -log(step_size)))
  vlim_upper <- -vlim_lower
  
  vval_lower <- f(vlim_lower, ...)
  if(any(vval_lower < 0))
    stop("There are negative values of f!")
  vflag_lower <- vval_lower > 0
  if(any(vflag_lower)) {
    if(vflag_lower[1])
      stop("f is already positive at maximum lower limit!")
    lower <- vlim_lower[c(vflag_lower[-1], TRUE)][1]
  }
  else {
    warning("No positive f value for lower limits!")
    lower <- vlim_lower[length(vlim_lower)]
  }
    
  vval_upper <- f(vlim_upper, ...)
  if(any(vval_upper < 0))
    stop("There are negative values of f!")
  vflag_upper <- vval_upper > 0
  if(any(vflag_upper)) {
    if(vflag_upper[1])
      stop("f is already positive at maximum upper limit!")
    upper <- vlim_upper[c(vflag_upper[-1], TRUE)][1]
  }
  else {
    warning("No positive f value for upper limits!")
    upper <- vlim_upper[length(vlim_upper)]
  }
  
  return(c(lower, upper))
}