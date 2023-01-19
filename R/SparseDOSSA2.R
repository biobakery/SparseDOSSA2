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
#' @param spike_metadata for metadata spike-in configurations. Must be one of two things:
#' a) , 
#' \itemize{
#' \item a character string of \code{"none"}, \code{"both"} \code{"abundance"}, 
#' or \code{"prevalence"}, indicating whether or not
#' association with metadata will be spiked in. For the spiked-in case, it
#' indicates if features' abundance/prevalence/both characteristics will be associated
#' with metadata (also see explanations for \code{metadata_effect_size} and 
#' \code{perc_feature_spiked_metadata})
#' \item a data.frame for detailed spike-in configurations. This is the more 
#' advanced approach, where detailed specification for metadata-microbial 
#' feature associations are provided. Note: if \code{spike_metadata} is provided
#' as a data.frame, then \code{metadata_matrix} must be provided as well 
#' (cannot be generated automatically). In this case, \code{spike_metadata} 
#' must have exactly four columns: \code{metadata_datum}, \code{feature_spiked}, 
#' \code{associated_property}, and \code{effect_size}. Each row of the data.frame 
#' configures one specific metadata-microbe association. Specifically:
#' \itemize{
#' \item \code{metadata_datum} (integer) indicates the column number for the metadata
#' variable to be associated with the microbe
#' \item \code{feature_spiked} (character) indicates the microbe name to be associated
#' with the metadata variable
#' \item \code{associated_property} (character, either \code{"abundance"} or 
#' \code{"prevalence"}), indicating the property of the microbe to be modified.
#' If you want the microbe to be associated with the metadata variable
#' in both properties, include two rows in \code{spike_metadata}, one for 
#' abundance and one for prevalence
#' \item \code{effect_size} (numeric) indicating the strength of the association. 
#' This corresponds to log fold change in non-zero abundance for
#' \code{"abundance"} spike-in, and log odds ratio for \code{"prevalence"} 
#' spike-in
#' }
#' }
#' @param metadata_effect_size (for when \code{spike_metadata} is \code{"abundance"},
#' \code{"prevalence"}, or \code{"both"}) effect size of the spiked-in associations. This is
#' non-zero log fold change for abundance spike-in, and log odds ratio for prevalence spike-in
#' @param perc_feature_spiked_metadata (for when \code{spike_metadata} is \code{"abundance"},
#' \code{"prevalence"}, or \code{"both"})
#' percentage of features to be associated with metadata 
#' @param metadata_matrix the user can provide a metadata matrix to use for spiking-in
#' of feature abundances. If using default (\code{NULL}) two variables will be generated: 
#' one continous, and a binary one of balanced cases and controls. Note: if
#' \code{spike_metadata} is provided as a data.frame, then the user must provide
#' \code{metadata_matrix} too
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
#' \code{spike_metadata} for the original \code{spike_metadata} parameter provided
#' by the user, \code{metadata_matrix} for the
#' metadata (either provided by the user or internally generated), and 
#' \code{feature_metadata_spike_df}
#' for detailed specification of which metadata variables were used to spike-in associations
#' with which features, in what properties at which effect sizes. This is the 
#' same as \code{spike_metadata} if the latter was provided as a data.frame.
#' }
#' }
#' @export
#' @author Siyuan Ma, \email{syma.research@gmail.com}
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
                         spike_metadata = "none",
                         metadata_effect_size = 1,
                         perc_feature_spiked_metadata = 0.05,
                         metadata_matrix = NULL,
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
    # check that Omega and Sigma should agree
    if(max(abs(Omega %*% template$EM_fit$fit$Sigma - 
               diag(rep(1, length(template$EM_fit$fit$pi0))))) > 1e-10)
       stop("Omega shoud be the inverse of Sigma!")
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
  
  # generate spiked-in association with metadata
  if(!(is.character(spike_metadata) | is.data.frame(spike_metadata)))
    stop("spike_metadata must be either character or a data.frame!")
  if(is.character(spike_metadata)) {
    spike_metadata <- match.arg(spike_metadata, 
                                choices = c("none", "abundance", 
                                            "prevalence", "both"))
  }
  if(identical(spike_metadata, "none")) {
    if(verbose)
      message("spike_metadata is \"none\", ",
              "no metadata association will be simulated...")
    
    mat_spiked_metadata <- mat_null
    feature_metadata_spike_df <- NULL
  } else {
    if(verbose)
      message("Spiking in metadata association...")
    
    # metadata_matrix
    if(is.null(metadata_matrix)) {
      if(is.data.frame(spike_metadata))
        stop("spike_metadata is provided as a data.frame. User must specify ",
             "metadata_matrix as well!")
      if(verbose)
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
    
    # metadata_effect_size
    if(length(metadata_effect_size) != 1 & 
       length(metadata_effect_size) != n_metadata)
      stop("Length of metadata_effect_size can only be either 1 or number of ",
           "columns of metadata_matrix!")
    if(length(metadata_effect_size) == 1)
      metadata_effect_size <- rep(metadata_effect_size, n_metadata)
    
    # feature_metadata_spike_df
    if(is.data.frame(spike_metadata)) {
      if(verbose) {
        message("spike_metadata is provided as a data.frame; ",
                "will use for simulating metadata association ",
                "(metadata_effect_size and perc_feature_spiked_metadata will ",
                "be ignored)...")
      }
      
      # check format
      if(!all(c("metadata_datum",
                "feature_spiked",
                "associated_property",
                "effect_size") %in%
              colnames(spike_metadata)))
        stop("spike_metadata does not follow the correct format! ",
             "Must have the following columns: metadata_datum, ",
             "feature_spiked, associated_property, and effect_size.")
      if(!all(spike_metadata$feature_spiked %in% 
              features))
        stop("feature_spiked in spike_metadata must provide the ",
             "spiked feature names!")
      if(!all(spike_metadata$metadata_datum %in% 
              seq(1, n_metadata)))
        stop("metadata_datum in spike_metadata must provide the ",
             "associated metadata column number!")
      if(!all(spike_metadata$associated_property %in% 
              c("prevalence", "abundance")))
        stop("associated_property in spike_metadata must be ",
             "either \"prevalence\" or \"abundance\"!")
      
      feature_metadata_spike_df <- spike_metadata
    } else {
      if(verbose)
        message("spike_metadata is specified as ", spike_metadata, "; ",
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
