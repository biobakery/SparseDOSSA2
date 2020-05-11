pi0_1 <- pi0_2 <- 0.8
rho <- 0.4
mat_rho <- matrix(c(1, rho, rho, 1), 2, 2)
a_samples <- rcopulasso(n = 100000, pi0 = c(pi0_1, pi0_2),
                        mu = c(0, 0), sigma = c(1, 1),
                        Omega = solve(mat_rho))
marginals <- get_marginals(a_samples)

Rho2(rho_p = rho, pi0_1 = pi0_1, pi0_2 = pi0_2,
     glim_1 = marginals$glim[1], glim_2 = marginals$glim[2],
     g0_1 = marginals$g0[1], g0_2 = marginals$g0[2],
     sigmaMod_1 = marginals$sigmaMod[1], sigmaMod_2 = marginals$sigmaMod[2])

data_g <- vapply(seq_len(ncol(a_samples)),
                 function(i_feature)
                   a_to_g(a = a_samples[, i_feature],
                          pi0 = marginals$pi0[i_feature],
                          mu = marginals$mu[i_feature],
                          sigma = marginals$sigma[i_feature]),
                 rep(0.0, nrow(a_samples)))
get_s(cor_data = cor(data_g),
      pi0 = marginals$pi0,
      glim = marginals$glim,
      g0 = marginals$g0,
      sigmaMod = marginals$sigmaMod)
