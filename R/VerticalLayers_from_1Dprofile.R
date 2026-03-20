#' Main ForestLayers function, use with 1D data
#'
#' @description Automatic layers delineation and parameterization from a 1D vertical profile.
#' Converts a 1D vertical profile (2 columns data.frame) into a list of the parameterized layers with their parameters and Weibull estimations.
#' The results are, for each layer, a data.frame with the initial VegDens and estimated VegDens, and a list of structural parameters.
#'
#' @param profile A 2 columns dataframe describing the VegDens vertical distribution (Z, VegDens). Can be the intermediate result of the VerticalLayers_from_3Dgrid() function.
#' @param diff_threshold numeric. The difference threshold (between 0 and 1) corresponding to a percentage of the segment’s maximum VegDens. Default is 0.25. See future paper for more details.
#' @param diff_max_threshold numeric. The difference threshold (between 0 and 1) corresponding to a percentage of the profile’s maximum VegDens. Default is 0.25. See future paper for more details.
#' @param col_name character. Name of the value of vegetation density (e.g. PAD, biomass, voxel count). Used for visualization purpose. Default is "VegDens".
#' @param print logical. Choose TRUE to visualize intermediate graphs created during the function.
#' @param debug logical. For development only.
#' @param export logical. Choose TRUE and specify a directory (output_dir) to export the 1D profile.
#' @param output_dir character. Output directory to export the intermediate results of the function.
#' @param plot_name character. Name of your data, for proper export procedure.
#' @return Returns, for each layer, a 3 columns data.frame with the height, the initial VegDens and the estimated VegDens, and a list of structural parameters. Returns also a visual representation of the full parameterized vertical profile.
#' @export
#'
#' @examples
#' \dontrun{
#' # Open your 2 columns data.frame and check that you have the correct columns (height, VegDens)
#' profile <- readxl::read_xlsx(my_file, sheet = "Sheet1")
#' profile <- as.data.frame(profile)
#'
#' # Define ForestLayers parameters
#' diff_threshold = 0.25
#' diff_max_threshold = 0.25
#' col_name <- "VegDens (units)"
#' output_dir <- "C:/my/output/directory"
#' plot_name <- "my_study_area"
#'
#' my_layers <- VerticalLayers_from_1Dprofile(
#' profile,
#' diff_threshold,
#' diff_max_threshold,
#' col_name,
#' print = FALSE,
#' debug = FALSE,
#' export = TRUE,
#' output_dir,
#' plot_name
#' )
#'
#' # Export results
#' results <- export_DVV(ForestLayersList = my_layers, output_dir, plot_name)
#' }
#'
#' @seealso Use VerticalLayers_from_3Dgrid() function if you have 3D data (grid, voxelized point cloud).
#' See export_DVV() function to export your results.
#' @keywords vertical profile, layer, Weibull
#' @author Amélie Juckler and other contributors
#'
#'

VerticalLayers_from_1Dprofile <- function(profile, diff_threshold = 0.25, diff_max_threshold = 0.25, col_name = "VegDens", print = FALSE, debug = FALSE, export = FALSE, output_dir = NULL, plot_name = NULL) {
  tryCatch({

    vertical_layers <- segments_PAD_profile(profile, diff_threshold, diff_max_threshold, col_name, print, debug, export, output_dir, plot_name)
    fitted_layers <- weibull_fit(segmentList = vertical_layers, col_name, print, debug)
    layersList = fitted_layers$weibullList
    final_layers <- final_limits(weibullList = layersList, col_name, print, debug, export, output_dir, plot_name)
    nrmse_final <- nrmse_global(finalList = final_layers)
    final_plot <- visu(finalList = final_layers, col_name, debug)
    finalResults <- get_final_results(finalList = final_layers)

    return(list(
      finalResults = finalResults,
      nrmse = nrmse_final[1],
      r2 = nrmse_final[2],
      final_plot = final_plot
    ))

  }, error = function(e) {
    stop(paste("Error in VerticalLayers_from_1Dprofile: ", e$message))
  })

}
