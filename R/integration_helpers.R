integrate2 <- function(f, 
                           lower, upper, 
                           rel_tol, abs_tol, max_eval, 
                           precBits,
                           ...) {
  if(lower == 0 & upper == 0)
    return(list(integral = 0,
                error = 0,
                neval = 0,
                returnCode = 1))
  
  ff <- function(x) f(x, ...)
  
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
    # linear spline for estimating integration
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
    
    # error estimation
    errors_spline <- estimate_errors(
      knots_diff, 
      c(Rmpfr::mpfr(0, precBits = precBits), 
        coefs_spline[2, ], 
        Rmpfr::mpfr(0, precBits = precBits)))
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
  return(list(integral = integral,
              error = error,
              neval = neval,
              returnCode = 0,
              debug = list(knots_spline = knots_spline,
                           vals_spline = vals_spline,
                           errors_spline = errors_spline)))
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

estimate_errors <- function(knots_diff, slopes) {
  if(length(knots_diff) != length(slopes) - 2)
    stop("Length of knots_diff and slopes should agree!")
  errors <- rep(NA, length = length(knots_diff))
  
  y3s <- knots_diff * slopes[-c(1, length(slopes))]
  errors <- 
    Rmpfr::sapplyMpfr(seq_along(knots_diff),
                      function(i_region)
                        estimate_one_error(
                          knots_diff[i_region],
                          y3s[i_region],
                          slopes[c(i_region, i_region + 1, i_region + 2)]
                        ))
  
  return(errors)
}

estimate_one_error <- function(x3, y3, slopes) {
  slope1 <- slopes[1]
  slope2 <- slopes[3]
  if((slope1 <= slopes[2] & slopes[2] <= slope2) |
     (slope1 >= slopes[2] & slopes[2] >= slope2)) {
    # function is convex/concave here
    
    # if slopes are the same
    if(slope1 == slope2)
      return(0)

    # if not calculate where lines meet
    x2 <- (slope2 * x3 - y3) / (slope2 - slope1)
    y2 <- x2 * slope1
    return(abs(x2 * y3 - x3 * y2) / 2)
  } else {
    return(abs(x3 * y3) / 2)
  }
}
