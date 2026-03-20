#' Main ForestLayers function, use with 3D data
#'
#' @description Automatic layers delineation and parameterization from a 3D matrix (X, Y, Z, VegDens).
#' Converts a 3D matrix (4 columns data.frame) into a list of the parameterized layers with their parameters and Weibull estimations.
#' The results are, for each layer, a data.frame with the initial VegDens and estimated VegDens, and a list of structural parameters.
#'
#' @param matrix 4 columns dataframe containing coordonates of each point with its VegDens value (X, Y, Z, VegDens).
#' @param dtm raster. Digital Terrain Model for the study area. The dtm will be opened with the function raster(). Necessary if normalize = TRUE.
#' @param vox_size numeric. Dimension of the voxels/grid (in meter). default = 0.3 m
#' @param plot_area numeric. Dimension of the matrix, the study area (in square meter). default = 1600 m2
#' @param normalize logical. Choose TRUE if you want to normalize the matrix and you have a raster dtm of the study area. Uses the lidR::height_above_ground() function.
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
#' # Open your 4 columns data.frame and check that you have the correct columns (x, y, z, VegDens)
#' matrix <- readxl::read_xlsx(my_file, sheet = "Sheet1")
#' matrix <- as.data.frame(profile)
#'
#' # Define ForestLayers parameters
#' dtm <- "c:/directory/of/my/dtm"
#' vox_size <- 0.3
#' plot_area <- 400
#' diff_threshold = 0.25
#' diff_max_threshold = 0.25
#' col_name <- "VegDens (units)"
#' output_dir <- "C:/my/output/directory"
#' plot_name <- "my_study_area"
#'
#' my_layers <- VerticalLayers_from_3Dgrid(
#' matrix,
#' dtm,
#' vox_size,
#' plot_area,
#' normalize = TRUE,
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
#' @seealso Use VerticalLayers_from_1Dprofile() function if you have 1D data (vertical profile).
#' See export_DVV() function to export your results.
#' @keywords layer, Weibull
#' @author Amélie Juckler and other contributors
#'
#'
VerticalLayers_from_3Dgrid <- function(matrix, dtm, vox_size = 0.3, plot_area = 400, normalize = T, diff_threshold = 0.25, diff_max_threshold = 0.25, col_name = "VegDens", print = F, debug = F, export = F, output_dir = NULL, plot_name = NULL) {
  tryCatch ({

    profile <- creates_PAD_profile(matrix, dtm, vox_size, plot_area, normalize, debug, export, output_dir, plot_name)
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
