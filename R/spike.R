spike_a_metadata <- function(null,
                             feature_param,
                             metadata,
                             spike_df) {
  if(ncol(null) != nrow(metadata))
    stop("Sample size of null abundance matrix and metadata do not agree!")
  if(nrow(spike_df) != 
     nrow(dplyr::distinct(spike_df, 
                          metadata_datum, 
                          feature_spiked, 
                          associated_property)))
    stop("feature-metadata spiking specification data frame has duplicate ",
         "metadata_datum, feature_spiked, and associated_property tuples!")
  spiked <- null
  for(feature_i in unique(spike_df$feature_spiked)) {
    spike_df_i <- subset(spike_df, feature_spiked == feature_i)
    spike_df_i_abundance <- subset(spike_df_i, 
                                   associated_property == "abundance")
    spike_df_i_prevalence <- subset(spike_df_i, 
                                    associated_property == "prevalence")
    spiked[feature_i, ] <- 
      spike_oneA_metadata(param = feature_param[feature_i, ],
                          metadata = metadata,
                          col_abundance = spike_df_i_abundance$metadata_datum,
                          effect_abundance = spike_df_i_abundance$effect_size,
                          col_prevalence = spike_df_i_prevalence$metadata_datum,
                          effect_prevalence = spike_df_i_prevalence$effect_size)
  }
  
  dimnames(spiked) <- dimnames(null)
  return(spiked)
}

spike_oneA_metadata <- function(param,
                                metadata,
                                col_abundance = c(),
                                effect_abundance = c(),
                                col_prevalence = c(),
                                effect_prevalence = c()) {
  effect_abundance_all <- rep(0, ncol(metadata))
  effect_prevalence_all <- rep(0, ncol(metadata))
  
  effect_abundance_all[col_abundance] <- effect_abundance
  effect_prevalence_all[col_prevalence] <- effect_prevalence
  
  pi0 <- expit(logit(param["pi0"]) - (metadata %*% effect_prevalence_all)[, 1])
  mu <- param["mu"] + (metadata %*% effect_abundance_all)[, 1]
  
  a <- rZILogN_one(n = nrow(metadata),
                   pi0 = pi0,
                   mu = mu,
                   sigma = param["sigma"])
  return(a)
}
