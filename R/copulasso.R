#' Title
#'
#' @param data 
#' @param lambda_list 
#' @param gamma_EBIC 
#' @param K_CV 
#'
#' @return
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
copulasso <- function(data, lambda_list,
                      gamma_EBIC = 0.5,
                      K_CV = 5) {

  df_eval <- data.frame(lambda = lambda_list,
                        df = NA_real_,
                        AIC = NA_real_,
                        BIC = NA_real_,
                        EBIC = NA_real_)
  s_data <- get_s(data = data)

  l_fits_evals <- future.apply::future_lapply(
    lambda_list,
    function(lambda) {
      fit_glasso <- glasso_wrapper(S = s_data,
                                   lambda = lambda,
                                   threshold = 1e-6)
      df <- (sum(fit_glasso != 0) - ncol(data)) / 2
      negLogLik <- negLogLik_mvn(S = s_data, Omega = fit_glasso)
      AIC <- negLogLik * nrow(s_data) + 2 * df
      BIC <- negLogLik * nrow(s_data) + log(nrow(data)) * df
      EBIC <- BIC + 4 * gamma_EBIC * log(ncol(data)) * df
      return(list(fit = fit_glasso,
                  evals = c(df, AIC, BIC, EBIC)))
    }
  )
  
  l_fits <- lapply(l_fits_evals, function(x) x$fit)
  df_eval[, c("df", "AIC", "BIC", "EBIC")] <-
    t(vapply(l_fits_evals,
             function(x) x$evals,
             rep(0.0, 4)))

  if(!is.null(K_CV)) {
    folds <- sample.int(n = K_CV, size = nrow(data), replace = TRUE)
    ## FIXME
    # parallelization not optimized
    # ideally should parallel over combinations of k and lambda
    # instead of just over k
    negLogLik_CV <- future.apply::future_lapply(
      sesq_len(K_CV),
      function(k) {
        data_train <- data[folds != k, ]
        s_train <- get_s(data = data_train)
        
        l_fit_glasso <- lapply(lambda_list,
                               function(lambda)
                                 glasso_wrapper(S = s_train,
                                                lambda = lambda,
                                                threshold = 1e-6)
        )
        
        data_test <- data[folds == k, ]
        s_test <- get_s(data = data_test)
        return(vapply(seq_along(lambda_list),
                      function(i)
                        negLogLik_mvn(S = s_test, Omega = l_fit_glasso[[i]]) *
                        nrow(s_test),
                      0.0))
      }
    ) %>% Reduce("+", .)
    
    df_eval$negLogLik_CV <- negLogLik_CV / nrow(data)
  }

  return(list(fits = l_fits, df_eval = df_eval))
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
