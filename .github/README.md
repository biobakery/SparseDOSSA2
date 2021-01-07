# "Simulating realistic microbial observations with SparseDOSSA2" 

Author Name: "Siyuan Ma"  
Affiliation: Harvard T.H. Chan School of Public Health.  
Broad Institute email: siyuan.ma@pennmedicine.upenn.edu 

# Introduction
SparseDOSSA2 an R package for fitting to and the simulation of realistic microbial abundance observations. It provides functionlaities for: a) generation of realistic synthetic microbial observations, b) spiking-in of associations with metadata variables for e.g. benchmarking or power analysis purposes, and c) fitting the SparseDOSSA 2 model to real-world microbial abundance observations that can be used for a). This vignette is intended to provide working examples for these functionalities.

```
library(SparseDOSSA2)
# tidyverse packages for utilities
library(magrittr)
library(dplyr)
library(ggplot2)
```

# Installation
SparseDOSSA2 is a Bioconductor package and can be installed via the following command.
```
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("SparseDOSSA2")
```
# Simulating realistic microbial observations with SparseDOSSA2
The most important functionality of SparseDOSSA2 is the simulation of realistic synthetic microbial observations. To this end, SparseDOSSA2 provides three pre-trained templates, "Stool", "Vaginal", and "IBD", targeting continuous, discrete, and diseased population structures.
```
Stool_simulation <- SparseDOSSA2(template = "Stool", 
                                 n_sample = 100, 
                                 n_feature = 100,
                                 verbose = TRUE)
Vaginal_simulation <- SparseDOSSA2(template = "Vaginal", 
                                   n_sample = 100, 
                                   n_feature = 100,
                                   verbose = TRUE)
```

# Fitting to microbiome datasets with SparseDOSSA2
SparseDOSSA2 provide two functions, fit_SparseDOSSA2 and fitCV_SparseDOSSA2, to fit the SparseDOSSA2 model to microbial count or relative abundance observations. For these functions, as input, SparseDOSSA2 requires a feature-by-sample table of microbial abundance observations. We provide with SparseDOSSA2 a minimal example of such a dataset: a five-by-five of the HMP1-II stool study.
```
data("Stool_subset", package = "SparseDOSSA2")
# columns are samples.
Stool_subset[1:2, 1, drop = FALSE]
```

## Fitting SparseDOSSA2 model with fit_SparseDOSSA2
fit_SparseDOSSA2 fits the SparseDOSSA2 model to estimate the model parameters: per-feature prevalence, mean and standard deviation of non-zero abundances, and feature-feature correlations. It also estimates joint distribution of these parameters and (if input is count) a read count distribution.
```
fitted <- fit_SparseDOSSA2(data = Stool_subset,
                           control = list(verbose = TRUE))
# fitted mean log non-zero abundance values of the first two features
fitted$EM_fit$fit$mu[1:2]
```

## Fitting SparseDOSSA2 model with fitCV_SparseDOSSA2
The user can additionally achieve optimal model fitting via fitCV_SparseDOSSA2. They can either provide a vector of tuning parameter values (lambdas) to control sparsity in the estimation of the correlation matrix parameter, or a grid will be selected automatically. fitCV_SparseDOSSA2 uses cross validation to select an "optimal" model fit across these tuning parameters via average testing log-likelihood. This is a computationally intensive procedure, and best-suited for users that would like accurate fitting to the input dataset, for best simulated new microbial observations on the same features as the input (i.e. not new features).
```
set.seed(1)
fitted_CV <- fitCV_SparseDOSSA2(data = Stool_subset,
                                         lambdas = c(0.1, 1),
                                         K = 2,
                                         control = list(verbose = TRUE))
# the average log likelihood of different tuning parameters
apply(fitted_CV$EM_fit$logLik_CV, 2, mean)
# The second lambda (1) had better performance in terms of log likelihood,
# and will be selected as the default fit.
```

# Parallelization controls with future
SparseDOSSA2 internally uses r BiocStyle::CRANpkg("future") to allow for parallel computation. The user can thus specify parallelization through future's interface. See the reference manual for future for more details. This is particularly suited if fitting SparseDOSSA2 in a high-performance computing environment/
```
## regular fitting 
# system.time(fitted_regular <- 
#               fit_SparseDOSSA2(data = Stool_subset,
#                                control = list(verbose = FALSE)))
## parallel fitting with future:
# future::plan(future::multisession())
# system.time(fitted_parallel <- 
#               fit_SparseDOSSA2(data = Stool_subset,
#                                control = list(verbose = FALSE)))

## For CV fitting, there are three components that can be paralleled, in order:
## different cross validation folds, different tuning parameter lambdas, 
## and different samples. It is usually most efficient to parallelize at the
## sample level:
# system.time(fitted_regular_CV <-
#               fitCV_SparseDOSSA2(data = Stool_subset,
#                                  lambdas = c(0.1, 1),
#                                  K = 2,
#                                  control = list(verbose = TRUE)))
# future::plan(future::sequential(), future::sequential(), future::multisession())
# system.time(fitted_parallel_CV <-
#               fitCV_SparseDOSSA2(data = Stool_subset,
#                                  lambdas = c(0.1, 1),
#                                  K = 2,
#                                  control = list(verbose = TRUE)))
```

# Sessioninfo
```
sessionInfo()

R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SparseDOSSA2_0.99.0 Rmpfr_0.8-2         gmp_0.6-1           igraph_1.2.6       
 [5] truncnorm_1.0-8     magrittr_2.0.1      future.apply_1.7.0  future_1.21.0      
 [9] huge_1.3.4.1        mvtnorm_1.1-1       ks_1.11.7           BiocCheck_1.22.0   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5          compiler_3.6.2      BiocManager_1.30.10 bitops_1.0-6       
 [5] tools_3.6.2         digest_0.6.27       mclust_5.4.7        jsonlite_1.7.2     
 [9] lattice_0.20-41     pkgconfig_2.0.3     Matrix_1.2-18       graph_1.64.0       
[13] curl_4.3            parallel_3.6.2      xfun_0.20           stringr_1.4.0      
[17] httr_1.4.2          knitr_1.30          globals_0.14.0      stats4_3.6.2       
[21] grid_3.6.2          getopt_1.20.3       optparse_1.6.6      Biobase_2.46.0     
[25] listenv_0.8.0       R6_2.5.0            parallelly_1.23.0   XML_3.99-0.3       
[29] RBGL_1.62.1         codetools_0.2-18    biocViews_1.54.0    BiocGenerics_0.32.0
[33] MASS_7.3-53         stringdist_0.9.6.3  RUnit_0.4.32        KernSmooth_2.23-18 
[37] stringi_1.5.3       RCurl_1.98-1.2     
```
