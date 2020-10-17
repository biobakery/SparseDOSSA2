integrate2 <- function(f, 
                       lower, upper, 
                       rel_tol, abs_tol, max_eval, 
                       precBits,
                       method, 
                       offset = FALSE, 
                       ...) {
  if(lower == 0 & upper == 0)
    return(list(integral = 0,
                error = 0,
                neval = 0,
                returnCode = 1))
  
  if(offset)
    val_offset <- f((lower + upper) / 2, ...)
  else
    val_offset <- 1
  
  ff <- function(x) f(x, ...) / val_offset 
  
  if (method %in% c("hcubature", "pcubature")) {
    fit_integrate <- cubature::cubintegrate(ff,
                                            lower = lower, upper = upper,
                                            relTol = rel_tol, absTol = abs_tol,
                                            method = method, maxEval = max_eval)
    fit_integrate$integral <- fit_integrate$integral * val_offset
    fit_integrate$error <- fit_integrate$error * val_offset
    return(fit_integrate)
  } else if (method == "spline") {
    
    neval <- 2
    knots_spline <- Rmpfr::mpfr(c(lower, upper),
                                precBits = precBits)
    vals_spline <- 
      Rmpfr::mpfr(
        ff(as.double(knots_spline)), 
        precBits = precBits)
    errors_spline <- Inf
    
    while(TRUE) {
      i_max_error <- which(errors_spline == max(errors_spline))[1]
      knots_spline <- c(knots_spline[seq(1, i_max_error)],
                        Rmpfr::mean(knots_spline[c(i_max_error, i_max_error + 1)]),
                        knots_spline[seq(i_max_error + 1, neval)])
      vals_spline <- c(vals_spline[seq(1, i_max_error)],
                       Rmpfr::mpfr(ff(as.double(knots_spline[i_max_error + 1])),
                                   precBits = precBits),
                       vals_spline[seq(i_max_error + 1, neval)])
      
      neval <- neval + 1
      knots_diff <-  knots_spline[-1] - knots_spline[-neval]
      errors_spline <- 
        abs(knots_diff *
              (vals_spline[-1] - vals_spline[-neval]))
      
      coefs_spline <- Rmpfr::mpfrArray(NA, precBits = precBits,
                                       dim = c(2, neval - 1))
      coefs_spline[2, ] <- 
        (vals_spline[-1] - vals_spline[-neval]) /
        knots_diff
      coefs_spline[1, ] <- 
        vals_spline[-neval] - knots_spline[-neval] * coefs_spline[2, ]
      
      integral <- sum(coefs_spline[1, ] * knots_spline[-1] +
                           coefs_spline[2, ] / 2 * knots_spline[-1]^2 -
                           coefs_spline[1, ] * knots_spline[-neval] - 
                           coefs_spline[2, ] / 2 * knots_spline[-neval]^2)
      
      error <- sum(errors_spline)
      
      if(neval >= max_eval)
        break
      if(integral < 0)
        stop("Negative integration values; something went wrong!")
      if(integral > 0)
        if(error / abs(integral) < rel_tol |
           error < abs_tol)
          break
    }
    return(list(integral = integral * val_offset,
                error = error * val_offset,
                neval = neval,
                returnCode = 0))
  }
}


get_intLimits <- function(x, pi0, mu, sigma, Omega, Sigma,
                          maxit = 10) {
  ind_nonzero <- x != 0
  logx <- log(x[ind_nonzero])
  mat_range <- cbind(mu[ind_nonzero] - logx - sqrt(2000)*sigma[ind_nonzero],
                     mu[ind_nonzero] - logx + sqrt(2000)*sigma[ind_nonzero])
  range <- c(max(c(mat_range[, 1], -745 - min(logx))), 
             min(mat_range[, 2]))
  
  if(range[1] >= range[2])
    return(c(0, 0)) ## FIXME??
  
  i_iter <- 1
  vlim <- seq(from = range[1], to = range[2], length.out = 3)
  vval <- sapply(vlim,
                 function(vv) 
                   dloga_asum(asum = exp(vv), x = x,
                              pi0 = pi0, mu = mu, sigma = sigma,
                              Omega = Omega, Sigma = Sigma))
  while(TRUE) {
    vflag <- vval > -745
    if(vflag[1] | rev(vflag)[1])
      stop("Positive dx values at integration limits!")
    if(sum(vflag) > 1)
      return(vlim[c(min(which(vflag)) - 1, 
                    max(which(vflag)) + 1)])
    
    if(i_iter + 1 > maxit)
      return(range(vlim)) ## FIXME??
    i_iter <- i_iter + 1
    
    if(all(vval == -Inf)) {
      vlim <- c(range[1],
                sort(runif(3, min = range[1], max = range[2])),
                range[2])
      vval <- c(-Inf,
                sapply(vlim,
                     function(vv)
                       dloga_asum(asum = exp(vv), x = x,
                                  pi0 = pi0, mu = mu, sigma = sigma,
                                  Omega = Omega, Sigma = Sigma)),
                -Inf)
    } else {
      ind_max <- order(-vval)[1]
      if(ind_max == 1) ind_max <- 2
      if(ind_max == length(vlim)) ind_max <- length(vlim) - 1
      vlim <- vlim[seq(ind_max - 1, 
                       ind_max + 1)]
      vval <- vval[seq(ind_max - 1, 
                       ind_max + 1)]
      
      vlim <- c(vlim[1],
                mean(vlim[seq(1, 2)]), 
                vlim[2], 
                mean(vlim[seq(2, 3)]),
                vlim[3])
      vval <- c(vval[1],
                dloga_asum(asum = exp(vlim[2]), x = x,
                           pi0 = pi0, mu = mu, sigma = sigma,
                           Omega = Omega, Sigma = Sigma),
                vval[2],
                dloga_asum(asum = exp(vlim[4]), x = x,
                           pi0 = pi0, mu = mu, sigma = sigma,
                           Omega = Omega, Sigma = Sigma),
                vval[3])
    }
  }
}

dloga_asum <- function(asum, x, pi0, mu, sigma, Omega, Sigma) {
  if(asum * min(x[x != 0]) <= 0)
    return(-Inf)
  if(asum > 0)
    return(dloga(a(x, asum), pi0 = pi0, mu = mu, sigma = sigma, 
                 Omega = Omega, Sigma = Sigma,
                 log.p = TRUE))
}

# get_intLimits2 <- function(x, pi0, mu, sigma, Omega, Sigma,
#                            control = list()) {
#   control <- do.call(control_integrate, control)
#   
#   limits <- get_intLimits(x = x, pi0 = pi0, mu = mu, sigma = sigma,
#                           Omega = Omega, Sigma = Sigma,
#                           maxit = control$maxit_limits)
#   
#   vlim <- seq(limits[1], limits[2], length.out = control$n_vals_limits)
#   vval <- sapply(vlim,
#                  function(vv)
#                    dloga(a = a(x, exp(vv)),
#                          pi0 = pi0, mu = mu, sigma = sigma,
#                          Omega = Omega, Sigma = Sigma))
#   vlim <- vlim[vval > -Inf]
#   vval <- vval[vval > -Inf]
#   
#   coef_quad <- lm(vval ~ vlim + I(vlim^2))$coef
#   mu_integrand <- - coef_quad[2] / coef_quad[3] / 2
#   sigma_integrand <- sqrt(- 1 / coef_quad[3] / 2)
#   
#   return(mu_integrand + 
#            qnorm(c(control$limit_tol / 2, 
#                    1 - control$limit_tol / 2)) * 
#            sigma_integrand)
#   
# }

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
