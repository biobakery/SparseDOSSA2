#' Create "count" template from Metaphlan output
#'
#' @param feature_abd Metaphlan feature abundance table
#' @param n_i Metaphlan estimate of read depth
#'
#' @return "Count" feature table that can be used as template for SparseDOSSA2
create_template_fromMetaphlan <- function(feature_abd, n_i) {
  if(ncol(feature_abd) != length(n_i))
    stop("Number of columns of feature_abd must be equal to the length of n_i!")
  mat_count <- round(t(t(feature_abd) * n_i))
  dimnames(mat_count) <- dimnames(feature_abd)
  return(mat_count)
}
