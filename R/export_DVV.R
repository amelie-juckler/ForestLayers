#' Export results
#'
#' @description Main function to export the final results of the package.
#'
#' @param ForestLayersList R list. Contains the return of VerticalLayers_from_1Dprofile() or VerticalLayers_from_3Dgrid() function.
#' @param output_dir character. Output directory to export the final results of the main function.
#' @param plot_name character. Name of your data, for proper export procedure.
#' @return Exports the results of layers delineation and parameterization into an Excel file + the visual representation of the vertical layers of your study area (png).
#' @export
#'
#' @seealso See VerticalLayers_from_1Dprofile() and VerticalLayers_from_3Dgrid() for utilization examples.
#' @keywords export
#' @author Amélie Juckler
#'
#'
export_DVV <- function(ForestLayersList, output_dir = NULL, plot_name = NULL) {

  cat("Exporting plot ", plot_name, "\n")
  cat("\n")

  layersList <- ForestLayersList$finalResults
  layersList <- lapply(layersList, as.data.frame)
  stats <- data.frame(NRMSE = ForestLayersList$nrmse, R2 = ForestLayersList$r2, stringsAsFactors = FALSE)
  layersList$quality_metrics <- stats
  write_xlsx(layersList, file.path(output_dir, paste0("ForestLayers_", plot_name, ".xlsx")))
  plot <- ForestLayersList$final_plot
  ggsave(file.path(output_dir, paste0("ForestLayers_visualization_", plot_name, ".png")), plot)

  cat("Plot ", plot_name, " exported", "\n")
  cat("\n")

}
