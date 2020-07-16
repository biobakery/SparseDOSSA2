#' Simulation of realistic microbial observations
#'
#' @param n_sample Number of samples to simulate
#' @param n_feature Number of features to simulate
#' @param template Template for simulation. This can be either one of the 
#' prescribed datasets or fitting results from SparseDOSSA2_fit
#' @param new_features Should SparseDOSSA2 simulate new microbial features, or 
#' the same as the template dataset?
#' @param spike_metadata Should SparseDOSSA2 simulate associations with 
#' metadata?
#' @param spike_strength Effect size of association with metadata
#' @param control 
#'
#' @return
#'
#' @export
#' @examples
SparseDOSSA2 <- function(template,
                                  n_sample,
                                  new_features = TRUE,
                                  n_feature,
                                  spike_metadata = c("none", 
                                                     "both",
                                                     "abundance",
                                                     "presence"),
                                  metadata_effect_size = 1,
                                  perc_feature_spiked_metadata = 0.05,
                                  metadata_matrix = NULL,
                                  feature_metadata_spike_df = NULL,
                                  median_read_depth = 50000,
                                  control = list()) {
  control <- do.call(control_generate, control)
  
  # generate per-feature params
  feature_param_template <- cbind("pi0" = template$EM_fit$fit$pi0,
                                  "mu" = template$EM_fit$fit$mu,
                                  "sigma" = template$EM_fit$fit$sigma)
  if(!new_features) {
    if(control$verbose) 
      message("new_features is FALSE, ",
              "adopt per-feature parameters from template ",
              "(n_feature will be ignored)...")
    n_feature <- nrow(feature_param_template)
    feature_param <- feature_param_template
    Omega <- template$EM_fit$fit$Omega
  } else {
    if(control$verbose) 
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
  if(control$verbose) 
    message("Generating null absolute abundance matrix...")
  mat_null <- generate_a(n = n_sample,
                         feature_param = feature_param,
                         Omega = Omega)
  
  # generate spiked-in associatio with metadata
  spike_metadata <- match.arg(spike_metadata)
  if(spike_metadata == "none") {
    if(control$verbose)
      message("spike_metadata is \"none\", ",
              "no metadata association will be simulated...")
    mat_spiked_metadata <- mat_null
  } else {
    if(control$verbose)
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
      if(control$verbose) {
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
      if(control$verbose)
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
    if(control$verbose)
      message("Generating feature abundances with spiked-in metadata ",
              "associations...")
    mat_spiked_metadata <- 
      spike_a_metadata(null = mat_null,
                       feature_param = feature_param,
                       metadata = metadata_matrix,
                       spike_df = feature_metadata_spike_df)
  }
  
  mat_rel <- apply(mat_spiked_metadata, 2, TSS)
  if(control$verbose) 
    message("Generating count matrix...")
  # generate read depth
  depth_new <- generate_depth(mu_depth = template$depth_fit["mu_depth"],
                              sigma_depth = template$depth_fit["sigma_depth"],
                              n = n_sample,
                              median_depth = median_read_depth)
  # generate read counts
  mat_count <- generate_count(rel = mat_rel,
                              depth = depth_new)
  
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
                                      feature_metadata_spike_df),
              control = control))
}

control_generate <- function(verbose = FALSE) {
  list(verbose = verbose)
}