#' Statistics computation
#'
#' @description NRMSE and R2 calculation for the full vertical profile.
#'
#' @param finalList R list. Contains the return of final_limits() function.
#' @return Final NRMSE and R2 of the Weibull estimations for the full vertical profile.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
nrmse_global <- function(finalList) {
  complet <- do.call(rbind, finalList)
  complet_no_sol <- complet[complet$z > 2, ] # delete values < 2m because it is often ground points
  global_nrmse <- (sqrt(mean((complet_no_sol$VegDens - complet_no_sol$VegDens_pred)^2, na.rm = TRUE)))/(max(complet_no_sol$VegDens) - min(complet_no_sol$VegDens))
  ss1 <- sum((complet_no_sol$VegDens - complet_no_sol$VegDens_pred)^2, na.rm = TRUE)
  ss2 <- sum((complet_no_sol$VegDens - mean(complet_no_sol$VegDens))^2, na.rm = TRUE)
  r2 <- 1 - (ss1/ss2)
  print(paste("Global NRMSE : ", global_nrmse))
  print(paste("Global R2 : ", r2))

  return(c(global_nrmse, r2))

}
