#' Title
#'
#' @param data 
#' @param lambda 
#' @param control 
#'
#' @return
#' @export
#'
#' @examples
fit_SparseDOSSA2 <- function(data, 
                             lambda = NULL,
                             control = list()) {
  # initialization, filtering
  if(any(data < 0))
    stop("data has negative values!")
  data <- t(data)
  control <- do.call(control_fit, control)
  if(control$verbose)
    message("Filtering for all-zero features/samples...")
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  if(!is.null(control$debug_dir))
    dir.create(control$debug_dir)
  
  # check if table is count and fit library size parameters if needed
  row_sums <- apply(data, 1, sum)
  if(all(row_sums > 1)) {
    if(control$verbose)
      message("Data appears to be count table. ",
              "Fitting library size distribution...")
    depth_fit <- fit_depth(depth = row_sums)
  } else {
    if(control$verbose)
      message("Data appears to be relative abundance table.")
    depth_fit <- c("mu_depth" = NA, "sigma_depth" = NA) ## FIXME
  } 
  
  # normalize data
  data <- t(apply(data, 1, TSS, correct = TRUE))
  
  # EM fitting
  if(control$verbose)
    message("Fitting EM algorithm...")
  EM_fit <- EM(data = data,
               lambda = lambda, 
               control = control)
  
  # fit per-feature parameter 3d Gaussian kernel density
  if(control$verbose)
    message("Fitting joint distribution of per-feature parameters...")
  F_fit <- fit_F(feature_param = cbind("pi0" = EM_fit$fit$pi0,
                                       "mu" = EM_fit$fit$mu,
                                       "sigma" = EM_fit$fit$sigma))
  
  return(list(EM_fit = EM_fit,
              F_fit = F_fit,
              depth_fit = depth_fit,
              l_filtering = l_filtering))
}

control_fit <- function(maxit = 100,
                        rel_tol = 1,
                        abs_tol = 1e-2,
                        control_numint = list(),
                        control_copulasso = list(),
                        verbose = FALSE,
                        debug_dir = NULL) {
  list(maxit = maxit,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       control_numint = control_numint,
       control_copulasso = control_copulasso,
       verbose = verbose,
       debug_dir = debug_dir)
}


#' Title
#'
#' @param data 
#' @param lambdas 
#' @param K 
#' @param control 
#'
#' @return
#' @export
#'
#' @examples
fitCV_SparseDOSSA2 <- function(data, 
                               lambdas,
                               K = 5,
                               control = list()) {
  # initialization, filtering
  if(any(data < 0))
    stop("data has negative values!")
  data <- t(data)
  control <- do.call(control_fit, control)
  if(control$verbose)
    message("Filtering for all-zero features/samples...")
  l_filtering <- filter_data(data)
  data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]
  if(!is.null(control$debug_dir))
    dir.create(control$debug_dir)
  
  # check if table is count and fit library size parameters if needed
  row_sums <- apply(data, 1, sum)
  if(all(row_sums > 1)) {
    if(control$verbose)
      message("Data appears to be count table. ",
              "Fitting library size distribution...")
    depth_fit <- fit_depth(depth = row_sums)
  } else {
    if(control$verbose)
      message("Data appears to be relative abundance table.")
    depth_fit <- c("mu_depth" = NA, "sigma_depth" = NA) ## FIXME
  } 
  
  # normalize data
  data <- t(apply(data, 1, TSS, correct = TRUE))
  
  # CV EM fitting
  if(control$verbose)
    message("Fitting EM algorithm with cross-validation...")
  EM_fit <- EM_CV(data = data,
                  lambdas = lambdas, 
                  K = K,
                  control = control)
  
  # fit per-feature parameter 3d Gaussian kernel density
  if(control$verbose)
    message("Fitting joint distribution of per-feature parameters...")
  F_fit <- fit_F(feature_param = cbind("pi0" = EM_fit$fit$pi0,
                                       "mu" = EM_fit$fit$mu,
                                       "sigma" = EM_fit$fit$sigma))
  
  return(list(EM_fit = EM_fit,
              F_fit = F_fit,
              depth_fit = depth_fit,
              l_filtering = l_filtering))
}