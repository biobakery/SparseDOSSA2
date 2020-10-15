get_mutilde <- function(data, eloga) {
  mutilde <- 
    vapply(seq_len(ncol(data)),
           function(j) {
             ind_j <- data[, j] > 0
             mean(eloga[ind_j] + 
                    log(data[ind_j, j]))
           },
           0.0)
  
  return(mutilde - mean(mutilde))
}

get_sigma2tilde <- function(data, eloga, eloga2) {
  mutilde <- 
    vapply(seq_len(ncol(data)),
           function(j) {
             ind_j <- data[, j] > 0
             mean(eloga[ind_j] + 
                    log(data[ind_j, j]))
           },
           0.0)
  
  sigma2tilde <- 
    vapply(seq_len(ncol(data)),
           function(j) {
             ind_j <- data[, j] > 0
             mean(eloga2[ind_j] + 
                    2 * eloga[ind_j] * log(data[ind_j, j]) +
                    log(data[ind_j, j])^2)
           },
           0.0)
  
  return(sigma2tilde - mutilde^2 + (mutilde - mean(mutilde))^2)
}

solve_mu <- function(mutilde, sigma2hat)
  return(mutilde - sigma2hat * sum(mutilde) / sum(sigma2hat))

solve_sigma2 <- function(mutilde, sigma2tilde, n_nonzero, 
                         control, maxiter = 1000) {
  sigma2hat_old <- 
    sigma2hat <- 
    sigma2tilde
  i_iter <- 0
  while(TRUE) {
    i_iter <- i_iter + 1
    for(j in seq_along(sigma2hat_old)) 
      # sigma2hat[j] <- solve_one_sigma2(
      #   sum(sigma2hat[-j] / n_nonzero[-j]),
      #   sum(mutilde)^2,
      #   sigma2tilde[j] - mutilde[j]^2,
      #   n_nonzero[j])
      sigma2hat[j] <- 
        sigma2tilde[j] - mutilde[j]^2
    # sigma2hat <- future.apply::future_vapply(
    #   seq_along(sigma2hat_old),
    #   function(j)
    #     solve_one_sigma2(
    #       sum(sigma2hat_old[-j]),
    #       sum(mutilde)^2,
    #       sigma2tilde[j] - mutilde[j]^2
    #     ),
    #   0.0
    # )
    diff_abs <- get_diff(sigma2hat, sigma2hat_old, 
                         denom_c = control$abs_tol / 10, method = "abs")
    diff_rel <- get_diff(sigma2hat, sigma2hat_old, 
                         denom_c = control$abs_tol / 10, method = "rel")
    
    if(max(diff_abs) < (control$abs_tol / 10) & 
       max(diff_rel) < (control$rel_tol / 10) |
       i_iter > maxiter) {
      break
    }
           
    sigma2hat_old <- sigma2hat
  }
  
  return(sigma2hat)
  
}

# solve_sigma22 <- function(mutilde, sigma2tilde, control, maxiter = 1000) {
#   sigma2hat_old <- 
#     sigma2hat <- 
#     sigma2tilde
#   i_iter <- 0
#   while(TRUE) {
#     i_iter <- i_iter + 1
#     for(j in sample.int(length(sigma2tilde))) 
#       sigma2hat[j] <- solve_one_sigma2(
#         sum(sigma2hat[-j]),
#         sum(mutilde)^2,
#         sigma2tilde[j] - mutilde[j]^2)
#     # sigma2hat <- future.apply::future_vapply(
#     #   seq_along(sigma2hat_old),
#     #   function(j)
#     #     solve_one_sigma2(
#     #       sum(sigma2hat_old[-j]),
#     #       sum(mutilde)^2,
#     #       sigma2tilde[j] - mutilde[j]^2
#     #     ),
#     #   0.0
#     # )
#     diff_abs <- get_diff(sigma2hat, sigma2hat_old, 
#                          denom_c = control$abs_tol / 10, method = "abs")
#     diff_rel <- get_diff(sigma2hat, sigma2hat_old, 
#                          denom_c = control$abs_tol / 10, method = "rel")
#     
#     if(max(diff_abs) < (control$abs_tol / 10) & 
#        max(diff_rel) < (control$rel_tol / 10) |
#        i_iter > maxiter) {
#       break
#     }
#     
#     sigma2hat_old <- sigma2hat
#   }
#   
#   return(sigma2hat)
#   
# }

solve_one_sigma2 <- function(a, b, c, n)
  uniroot(f = function(x) - n^2 * x * (x / n + a)^2 + x^2 * b + n^2 * c * (x / n + a)^2,
           interval = c(c, b + c))$root

