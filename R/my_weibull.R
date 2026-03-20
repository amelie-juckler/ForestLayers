#' Weibull function
#'
#' @param params List of Weibull parameters (4).
#' @param x Height values.
#' @return Returns the VegDens estimations with Weibull parameters for each x value.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
my_weibull <- function(params, x) {
  k <- params[1]
  l <- params[2]
  scale_factor <- params[3] # Otherwise need to normalize y between 0 and 1
  offset <- params[4] # Otherwise y values need to start at 0
  return(offset + scale_factor * (k / l) * ((x / l)^(k - 1)) * exp(-(x / l)^k)) # Weibull function
}
