test_that("spike in works", {
  # # no spike-in
  # sim <- SparseDOSSA2()
  # 
  # # spike-in default
  # sim <- SparseDOSSA2(spike_metadata = "both")
  # plot(sim$spike_metadata$metadata_matrix[, 1],
  #      sim$simulated_matrices$rel[60,])
  # boxplot(sim$simulated_matrices$rel[19, ] ~
  #           sim$spike_metadata$metadata_matrix[, 2])
  # 
  # # spike-in matrix
  # # error
  # spike_metadata <- sim$spike_metadata$feature_metadata_spike_df[c(1, 20), ]
  # sim <- SparseDOSSA2(spike_metadata = spike_metadata)
  # 
  # # error
  # metadata_matrix <- matrix(rnorm(100), ncol = 1)
  # sim <- SparseDOSSA2(spike_metadata = spike_metadata,
  #                     metadata_matrix = metadata_matrix)
  # 
  # # runs
  # metadata_matrix <- matrix(rnorm(300), ncol = 3)
  # sim <- SparseDOSSA2(spike_metadata = spike_metadata,
  #                     metadata_matrix = metadata_matrix)
})