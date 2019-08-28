#' correlation estimation with randomized option to break ties
#'
#' @param x numeric matrix x
#' @param y numeric matrix y (set to x if unprovided)
#' @param method method used for correlation estimation. Only "spearman" and "kendall" are supported.
#' @param random if randomization should be performed in order to break ties
#' @param R number of random iterations to perform
#'
#' @return estimated correlation matrix
cor2 <- function(x, y = x, method = "spearman",
                 random = FALSE, R = 30) {
  if(nrow(x) != nrow(y))
    stop("Number of rows between x and y must agree!")
  if(!(all(x >= 0) & all(y >= 0)))
    stop("Only suitable for feature abundance (i.e., non-negative) table!")

  method <- match.arg(method, c("spearman", "kendall"))

  # If not random just normal correlation
  if(!random) {
    return(cor(x, y, method = method))
  }

  # If randomize then replace zeros with randomized small values
  if(random) {
    ind_x <- x != 0
    ind_y <- y != 0
    minMeasure_x <- min(setdiff(abs(x), 0)) / 2
    minMeasure_y <- min(setdiff(abs(y), 0)) / 2
    x_fill <- x
    y_fill <- y
    l_cor <- lapply(1:R, function(r) {
      x_fill[!ind_x] <- runif(n = sum(!ind_x), min = -minMeasure_x, max=minMeasure_x)
      y_fill[!ind_y] <- runif(n = sum(!ind_y), min = -minMeasure_y, max=minMeasure_y)
      return(cor(x_fill, y_fill, method = method))
    })
    return(Reduce("+", l_cor)/R)
  }
}
