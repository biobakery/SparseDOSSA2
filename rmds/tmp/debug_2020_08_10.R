library(magrittr)
smar::sourceDir("R/")
# future::plan(future::multisession())

load("../data/physeqs/HMP1II_vaginal.RData")
x_samples_vaginal <- physeq_vaginal %>% 
  smar::otu_table2() %>% 
  t()
# subset to top 200 samples and 250 features
x_samples_vaginal <- 
  x_samples_vaginal[order(-apply(x_samples_vaginal, 1, sum))[seq(1, 200)], ]
x_samples_vaginal <- 
  x_samples_vaginal[, apply(x_samples_vaginal > 0, 2, sum) >= 2]
x_samples_vaginal <- 
  x_samples_vaginal[, apply(x_samples_vaginal, 1, TSS) %>% 
                      apply(1, mean) %>% 
                      {order(-.)[1:250]}]
# set to relative abundance
x_samples_vaginal <- 
  x_samples_vaginal %>% 
  apply(1, TSS) %>% 
  t()

data = x_samples_vaginal
lambdas = 1e-3
K = 10
control = list(abs_tol = 1e-2,
               rel_tol = 1e-2,
               maxit = 100,
               verbose = TRUE)

control <- do.call(control_fit, control)
l_filtering <- filter_data(data)
data <- data[l_filtering$ind_sample, l_filtering$ind_feature, drop = FALSE]


control <- do.call(control_fit, control)
load("../results/Stool_vaginal/1/K_0/lambda_1/debug_EM.RData")
l_fits_full <- list(list(fit = rev(l_debug$ll_params)[[1]],
                         lambda = lambdas[1],
                         convergence = list(converge = TRUE,
                                            converge_code = 0,
                                            n_iter = length(l_debug$ll_params)),
                         l_filtering = l_filtering))
load("../results/Stool_vaginal/1/CV_folds.RData")

l_logLik <- list()
for(k in seq_len(K)) {
  data_training <- data[CV_folds != k, ]
  data_testing <- data[CV_folds == k, ]
  
  load(paste0("../results/Stool_vaginal/1/K_", k, "/lambda_1/debug_EM.RData"))
  l_filtering <- filter_data(data_training)
  l_fits_k <- 
    list(list(fit = rev(l_debug$ll_params)[[1]],
              lambda = lambdas[1],
              convergence = list(converge = TRUE,
                                 converge_code = 0,
                                 n_iter = length(l_debug$ll_params)),
              l_filtering = l_filtering))
  
  # Fill in parameters estimates for features not present in training data
  for(i_lambda in seq_along(lambdas))
    l_fits_k[[i_lambda]]$fit <-
    fill_estimates_CV(l_fits_k[[i_lambda]]$fit,
                      l_fits_full[[i_lambda]]$fit,
                      l_fits_k[[i_lambda]]$l_filtering$ind_feature)
  
  # Calculate ll in testing data
  params <- l_fits_k[[1]]$fit
  l_logLik[[k]] <-future.apply::future_sapply(
    seq_len(nrow(data_testing)),
    function(i_sample) {
      logLik <-
        dx(x = data_testing[i_sample, , drop = TRUE],
           pi0 = params$pi0, mu = params$mu, sigma = params$sigma,
           Omega = params$Omega, Sigma = params$Sigma,
           control = control$control_numint,
           log.p = TRUE)
      logLik
    })
}

for(i_sample in seq_len(nrow(data_testing)))
  logLik <-
  dx(x = data_testing[i_sample, , drop = TRUE],
     pi0 = params$pi0, mu = params$mu, sigma = params$sigma,
     Omega = params$Omega, Sigma = params$Sigma,
     control = control$control_numint,
     log.p = TRUE)
