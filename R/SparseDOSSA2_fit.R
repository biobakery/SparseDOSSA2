#' Fit SparseDOSSA 2 model to a microbiome abundance dataset
#' 
#' \code{fit_SparseDOSSA2} fits the SparseDOSSA 2 model (zero-inflated log normal
#' marginals connected through Gaussian copula) to microbial abundances. It takes
#' as input a feature-by-sample microbial count or relative abundance table and
#' a penalization tuning parameter \code{lambda} to control the sparsity of 
#' feature-feature correlations. It then adopts a penalized expectation-maximization
#' algorithm to provide estimations of the model parameters.
#' 
#' @param data feature-by-sample matrix of abundances (proportions or
#' counts)
#' @param lambda positive penalization parameter for the sparsity of feature-feature
#' correlations. Default to maxmum value \code{1}, where features are assumed to be 
#' independent (no correlations, most sparse)
#' @param control a named list of additional control parameters. See help page for
#' \code{control_fit}
#'
#' @return a list, with the following components:
#' \describe{
#' \item{EM_fit}{
#' list of fitted parameters from the EM algorithm. 
#' }
#' \item{F_fit}{
#' fitted parameters for the joint distribution of per-feature prevalence, abundance,
#' and variability parameters (for simulating new features) 
#' }
#' \item{depth_fit}{fitted parameters for the read depth distribution. Only applicable
#' to count data.
#' }
#' \item{l_filtering}{list of quality control filtering for sample and features.
#' }
#' }

#' @export
#' @author Siyuan Ma, \email{siyuan.ma@pennmedicine.upenn.edu}
#' @examples
#' data("Stool_subset")
#' fitted <- fit_SparseDOSSA2(data = Stool_subset)
fit_SparseDOSSA2 <- function(data, 
                             lambda = 1,
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

#' Control options for fit_SparseDOSSA2 and fitCV_SparseDOSSA2
#'
#' @param maxit maximum number of EM iterations
#' @param rel_tol relative change threshold in the log likelihood
#' for algorithm convergence
#' @param abs_tol absolute change threshold in the log likelihood
#' for algorithm convergence
#' @param control_numint a named list of control parameters for the 
#' numerical integrations during the E step. See help page for
#' \code{control_numint}
#' @param verbose whether or not detailed running messages should be provided
#' @param debug_dir directory for intermediate output, such as the
#' EM expectations and parameter values and during each step of the
#' EM algorithm. Default to \code{NULL} in which case no such output
#' will be generated
#'
#' @return a list of the same names
#' @export
control_fit <- function(maxit = 100,
                        rel_tol = 1e-2,
                        abs_tol = 1e-2,
                        control_numint = list(),
                        verbose = FALSE,
                        debug_dir = NULL) {
  list(maxit = maxit,
       rel_tol = rel_tol,
       abs_tol = abs_tol,
       control_numint = control_numint,
       verbose = verbose,
       debug_dir = debug_dir)
}


#' Fit SparseDOSSA 2 model to a microbiome abundance dataset with cross validation
#' 
#' \code{fitCV_SparseDOSSA2} randomly partitions the data into fitting and testing
#' subsets. It fits the SparseDOSSA 2 model to the fitting sets and uses log likelihood
#' of the fitted parameters in the testing sets as the criteria for selection of 
#' tuning parameter lambda.
#' 
#' @param data feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param lambdas vector of positive penalization parameters for the sparsity of feature-feature
#' correlations. The function fits SparseDOSSA 2 models to each of the lambda values, and uses
#' cross validation likelihood to select the optimal one. If not provided this will be chosen
#' automatically.
#' @param control a named list of additional control parameters. See help page for
#' \code{control_fit}.
#'
#' @return a list, with the following components:
#' \describe{
#' \item{EM_fit}{
#' list of fitted parameters from the EM algorithm, with additional cross validation likelihood. 
#' }
#' \item{F_fit}{
#' fitted parameters for the joint distribution of per-feature prevalence, abundance,
#' and variability parameters (for simulating new features) 
#' }
#' \item{depth_fit}{fitted parameters for the read depth distribution. Only applicable
#' to count data.
#' }
#' \item{l_filtering}{list of quality control filtering for sample and features.
#' }
#' }

#' @export
#' @author Siyuan Ma, \email{siyuan.ma@pennmedicine.upenn.edu}
#' @examples
#' data("Stool_subset")
#' fitted <- fitCV_SparseDOSSA(data = Stool_subset,
#'                             lambdas = c(0.1, 1),
#'                             K = 5)
#' 
fitCV_SparseDOSSA2 <- function(data, 
                               lambdas = 10^seq(-2, 0, length.out = 5),
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