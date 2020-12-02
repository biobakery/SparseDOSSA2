#' Simulate synthetic microbial abundance observations with SparseDOSSA2 
#'
#' \code{SparseDOSSA2} generates synthetic microbial abundance observations
#' from either pre-trained template, or user-provided fitted results from
#' \code{fit_SparseDOSSA2} or \code{fitCV_SparseDOSSA2}. Additional options
#' are available for simulating associations between microbial features
#' and metadata variables.
#' 
#' @param template can be 1) a character string (\code{"Stool"}, \code{"Vaginal"},
#' or \code{"IBD"}) indicating one of the pre-trained templates in SparseDOSSA2,
#' or 2) user-provided, fitted results. In the latter case this should be an output
#' from \code{fit_SparseDOSSA2} or \code{fitCV_SparseDOSSA2}.
#' @param n_sample number of samples to simulate
#' @param new_features \code{TRUE}/\code{FALSE} indicator for whether or not new
#' features should be simulated. If \code{FALSE} then the same set of features 
#' in \code{template} will be simulated.
#' @param n_feature number of features to simulate. Only relevant when 
#' \code{new_features} is \code{TRUE}
#' @param spike_metadata character string of \code{"none"}, \code{"both"}, 
#' \code{"abundance"}, or \code{"prevalence"}, indicating whether or not
#' association with metadata will be spiked in. For the spiked-in case, it also
#' indicates if features' abundance/prevalence/both characteristics will be associated
#' with metadata
#' @param metadata_effect_size effect size of the spiked-in associations. This is
#' log fold change for abundance spike-in, and log odds ratio for prevalence spike-in
#' @param perc_feature_spiked_metadata percentage of features to be associated with metadata 
#' @param metadata_matrix the user can provide a metadata matrix to use for spiking-in
#' of feature abundances. If using default (\code{NULL}) a single variable of balanced
#' cases and controls will be generated
#' @param feature_metadata_spike_df the user can provide a data frame to exactly prescribe
#' feature-metadata associations. This data frame should have exactly the following columns:
#' a) \code{"metadata_datum"} for which metadata to be associated with
#' b) \code{"feature_spiked} which features, in
#' c) \code{"associated_property"} which property (\code{"abundance"}, \code{"prevalence"}, or \code{"both"}), at
#' d) \code{"effect_size"} what effect size. Default to \code{NULL} in which case this
#' data frame will be automatically generated for the user, based on \code{spike_metadata}, 
#' \code{metadata_effect_size}, \code{perc_feature_spiked_metadata}, and \code{metadata_matrix} (if provided).
#' @param median_read_depth targeted median per-sample read depth 
#' @param verbose whether detailed information should be printed
#'
#' @return a list with the following component:
#' \describe{
#' \item{simulated_data}{
#' feature by sample matrix of simulated microbial count observations
#' }
#' \item{simulated_matrices}{
#' list of all simulated data matrices, including that of null (i.e. not spiked-in) absolute
#' abundances, spiked-in absolute abundances, and normalized relative abundances
#' }
#' \item{params}{
#' parameters used for simulation. These are provided in \code{template}.
#' }
#' \item{spike_metadata}{
#' list of variables provided or generated for metadata spike-in. This include 
#' \code{spike_metadata} for the characteristic specified (\code{"abundance"},  
#' \code{"prevalence"}, or \code{"both}), \code{metadata_matrix} for the
#' metadata (either provided by the user or internally generated), and 
#' \code{feature_metadata_spike_df} (either provided by the user or internally generated) 
#' for detailed specification of which metadata variables were used to spike-in associations
#' with which features, in what properties at which effect sizes.
#' }
#' }
#' @export
#' @author Siyuan Ma, \email{siyuan.ma@pennmedicine.upenn.edu}
#'
#' @examples
#' ## Using one of the pre-trained SparseDOSSA2 templates:
#' sim <- SparseDOSSA2(template = "stool", n_sample = 200, new_features = FALSE)
#' ## Using user-provided trained SparseDOSSA2 model:
#' data("Stool_subset")
#' fitted <- fit_SparseDOSSA(data = Stool_subset)
#' sim <- SparseDOSSA2(template = fitted, n_sample = 200, new_features = FALSE)
SparseDOSSA2 <- function(template = "Stool",
                         n_sample = 100,
                         new_features = TRUE,
                         n_feature = 100,
                         spike_metadata = c("none", 
                                            "both",
                                            "abundance",
                                            "prevalence"),
                         metadata_effect_size = 1,
                         perc_feature_spiked_metadata = 0.05,
                         metadata_matrix = NULL,
                         feature_metadata_spike_df = NULL,
                         median_read_depth = 50000,
                         verbose = TRUE) {
  if(is.character(template)) {
    if(!template %in% c("Stool", "Vaginal", "IBD"))
      stop("Pre-trained template must be one of \"Stool\", \"Vaginal\", or \"IBD\"!")
    template <- get(template)
  }
  
  # generate per-feature params
  feature_param_template <- cbind("pi0" = template$EM_fit$fit$pi0,
                                  "mu" = template$EM_fit$fit$mu,
                                  "sigma" = template$EM_fit$fit$sigma)
  if(!new_features) {
    if(verbose) 
      message("new_features is FALSE, ",
              "adopt per-feature parameters from template ",
              "(n_feature will be ignored)...")
    n_feature <- nrow(feature_param_template)
    feature_param <- feature_param_template
    Omega <- template$EM_fit$fit$Omega
  } else {
    if(verbose) 
      message("new_features is TRUE, ",
              "generating new per-feature parameters, based on template...")
    feature_param <- 
      generate_featureParam(F_fit = template$F_fit,
                            param_original = feature_param_template,
                            n_feature = n_feature)
    Omega <- diag(rep(1, nrow(feature_param)))
  }
  features <- 
    rownames(Omega) <- 
    colnames(Omega) <- 
    rownames(feature_param)
  
  # generate null absolute abundance matrix
  if(verbose) 
    message("Generating null absolute abundance matrix...")
  mat_null <- generate_a(n = n_sample,
                         feature_param = feature_param,
                         Omega = Omega)
  
  # generate spiked-in associatio with metadata
  spike_metadata <- match.arg(spike_metadata)
  if(spike_metadata == "none") {
    if(verbose)
      message("spike_metadata is \"none\", ",
              "no metadata association will be simulated...")
    mat_spiked_metadata <- mat_null
  } else {
    if(verbose)
      message("Spiking in metadata association...")
    if(is.null(metadata_matrix)) {
      message("metadata_matrix is not provided; ",
              "simulating default metadata_matrix...")
      metadata_matrix <- cbind(rnorm(n = n_sample),
                               rbinom(n = n_sample,
                                      size = 1,
                                      prob = 0.5))
      rownames(metadata_matrix) <- colnames(mat_null)
    } else {
      if(!is.matrix(metadata_matrix))
        stop("metadata_matrix must be a matrix ",
             "(model matrix where categorical variables are dummified)!")
      if(nrow(metadata_matrix) != n_sample)
        stop("n_sample does not agree with number of samples in ",
             "metadata_matrix!")
      if(!is.null(rownames(metadata_matrix)))
        colnames(mat_null) <- rownames(metadata_matrix)
    }
    n_metadata <- ncol(metadata_matrix)
    if(length(metadata_effect_size) != 1 & 
       length(metadata_effect_size) != n_metadata)
      stop("Length of metadata_effect_size can only be either 1 or number of ",
           "columns of metadata_matrix!")
    if(length(metadata_effect_size) == 1)
      metadata_effect_size <- rep(metadata_effect_size, n_metadata)
    if(!is.null(feature_metadata_spike_df)) {
      if(verbose) {
        message("feature_metadata_spike_df is provided; ",
                "will use for simulating metadata association ",
                "(metadata_effect_size and perc_feature_spiked_metadata will ",
                "be ignored)...")
      }
      if(!is.data.frame(feature_metadata_spike_df))
        stop("feature_metadata_spike_df must be a data frame!")
      if(!all(c("metadata_datum",
                "feature_spiked",
                "associated_property",
                "effect_size") %in%
              colnames(feature_metadata_spike_df)))
        stop("feature_metadata_spike_df does not follow the correct format! ",
             "Must have the following columns: metadata_datum, ",
             "feature_spiked, associated_property, and effect_size.")
      if(!all(feature_metadata_spike_df$feature_spiked %in% 
              features))
        stop("feature_spiked in feature_metadata_spike_df must provide the ",
             "spiked feature names!")
    } else {
      if(verbose)
        message("feature_metadata_spike_df is not provided; ",
                "generating default metadata association...")
      feature_metadata_spike_df <- 
        generate_feature_metadata_spike_df(
          features = features,
          perc_feature_spiked_metadata = perc_feature_spiked_metadata,
          n_metadata = n_metadata,
          effect_size = metadata_effect_size,
          spike_metadata = spike_metadata) 
    }
    if(verbose)
      message("Generating feature abundances with spiked-in metadata ",
              "associations...")
    mat_spiked_metadata <- 
      spike_a_metadata(null = mat_null,
                       feature_param = feature_param,
                       metadata = metadata_matrix,
                       spike_df = feature_metadata_spike_df)
  }
  
  mat_rel <- apply(mat_spiked_metadata, 2, TSS)
  if(any(is.na(template$depth_fit))) {
    mat_count <- mat_rel
  } else {
    if(verbose) 
      message("Generating count matrix...")
    # generate read depth
    depth_new <- generate_depth(mu_depth = template$depth_fit["mu_depth"],
                                sigma_depth = template$depth_fit["sigma_depth"],
                                n = n_sample,
                                median_depth = median_read_depth)
    # generate read counts
    mat_count <- generate_count(rel = mat_rel,
                                depth = depth_new)
  }
  
  return(list(simulated_data = mat_count,
              simulated_matrices = list(rel = mat_rel,
                                        a_spiked = mat_spiked_metadata,
                                        a_null = mat_null),
              params = list(feature_param = feature_param,
                            Omega = Omega),
              template = template,
              spike_metadata = list(spike_metadata = spike_metadata,
                                    metadata_matrix = metadata_matrix,
                                    feature_metadata_spike_df = 
                                      feature_metadata_spike_df)))
}