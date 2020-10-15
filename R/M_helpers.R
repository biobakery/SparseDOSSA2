get_muhat <- function(data, eloga) {
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

get_sigmahat <- function(data, eloga, eloga2) {
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
  
  return(sqrt(sigma2tilde - mutilde^2))
}