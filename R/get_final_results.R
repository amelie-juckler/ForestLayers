#' Formating results
#'
#' @description Function that gives the results of the main function in their final format.
#' Dividing the data.frame into a data.frame with the initial VegDens and estimated VegDens, and a list of structure parameters.
#'
#' @param finalList R list. Contains the return of final_limits() function.
#' @return Returns, for each layer, a 3 columns data.frame with the height, the initial VegDens and the estimated VegDens, and a list of structure parameters.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
get_final_results <- function(finalList) {

  tryCatch({

    results <- list()

    for (i in 1:(length(finalList))) {

      segment <- finalList[[i]]

      segment_final <- data.frame(
        z = segment$z,
        VegDens = segment$VegDens,
        VegDens_pred = segment$VegDens_pred
      )

      parameters_final <- list(
        K = segment$K[1],
        L = segment$L[1],
        Scale_Factor = segment$Scale_Factor[1],
        Offset = segment$Offset[1],
        VegDens_max_init = segment$VegDens_max_init[1],
        VegDens_max_weibull = segment$VegDens_max_weibull[1],
        peak_position_init = segment$peak_position_init[1],
        peak_position_weibull = segment$peak_position_weibull[1],
        NRMSE = segment$NRMSE[1],
        total_VegDens_init = segment$total_VegDens_init[1],
        total_VegDens_weibull = segment$total_VegDens_weibull[1],
        VegDens_mean_init = segment$VegDens_mean_init[1],
        VegDens_mean_weibull = segment$VegDens_mean_weibull[1],
        VegDens_med_init = segment$VegDens_med_init[1],
        VegDens_med_weibull = segment$VegDens_med_weibull[1],
        VegDens_var_init = segment$VegDens_var_init[1],
        VegDens_var_weibull = segment$VegDens_var_weibull[1],
        quantile_05 = segment$quantile_05[1],
        quantile_25 = segment$quantile_25[1],
        quantile_75 = segment$quantile_75[1],
        quantile_95 = segment$quantile_95[1],
        skewness = segment$skewness[1]
      )

      results[[paste0("Layer ", i)]] <- segment_final

      results[[paste0("Weibull parameters ", i)]] <- parameters_final

    }

  }, error = function(e){
    message("Warning detected during the execution of get_final_results function", e$message)
    results <- list()
  })

  return(results)

}
