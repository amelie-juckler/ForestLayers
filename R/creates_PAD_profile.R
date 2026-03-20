#' Conversion to 1D profile
#'
#' @description Converts a 3D matrix into a 1D vertical profile.
#' 3D matrix is a 4 columns dataframe describing the VegDens distribution (X, Y, Z, VegDens).
#'
#' @param matrix A 4 columns dataframe containing coordonates of each point with its VegDens value (X, Y, Z, VegDens).
#' @param dtm raster. Digital Terrain Model for the study area. Necessary if normalize = TRUE.
#' @param vox_size numeric. Dimension of the voxels/grid (in meter). default = 0.3 m
#' @param plot_area numeric. Dimension of the matrix, the study area (in square meter). default = 1600 m2
#' @param normalize logical. Choose TRUE if you want to normalize the matrix and you have a raster dtm of the study area. Uses the height_above_ground() function of lidR package.
#' @param debug logical. For development only.
#' @param export logical. Choose TRUE and specify a directory (output_dir) to export the 1D profile.
#' @param output_dir character. Output directory to export the intermediate results of the function.
#' @param plot_name character. Name of your data, for proper export procedure.
#' @return Returns a 2 columns dataframe (1D profile) describing the VegDens vertical distribution (z, VegDens).
#' @keywords internal
#' @author Amélie Juckler
#'
#'
creates_PAD_profile <- function(matrix, dtm = NULL, vox_size = 0.3, plot_area = 1600, normalize = F, debug = F, export = F, output_dir = NULL, plot_name = NULL) {

  tryCatch({

    profile <- data.frame()

    # Check column names
    stopifnot("X" %in% colnames(matrix) & "Y" %in% colnames(matrix))
    name <- names(matrix)[4]
    grid_df <- dplyr::rename(matrix, x = X, y = Y, Z = Z, PAD = name)

    if (debug) {
      print("dataframe has the right column names")
    }

    # Normalization
    if (normalize) {

      vox = LAS(data.frame(grid_df, x = grid_df$x, y = grid_df$y, z = grid_df$Z))
      vox@data <- vox@data[, c("X", "Y", "Z", "PAD")]
      mnt = raster(dtm)

      if (debug) {
        print("dtm found")
      }

      vox <- height_above_ground(vox, mnt)
      vox <- filter_poi(vox, vox@data$PAD > 0)
      vox@data$Classification[vox@data$hag < 0.5] <- as.integer(0) # Ground points (class 0)
      vox@data$Classification[vox@data$hag >= 0.5] <- as.integer(2) # Vegetation points (class 2)

      if (debug) {

        x <- plot(vox, color = "Classification", bg = "black")
        add_dtm3d(x, mnt)
        print("affichage 3d - ok")

      }

      vox1 <- filter_poi(vox, vox@data$Classification != 0) # Remove ground points

      grid_df <- as.data.frame(vox1@data)
      grid_df <- grid_df %>% rename(x = X, y = Y, z0 = Z, PAD = PAD, z = hag)
      grid_df <- grid_df %>% mutate(Z = floor(z/vox_size)*vox_size + (vox_size/2)) # revoxelization

    }

    # Creation of the 1D profile
    sum_VegDens <- aggregate(PAD ~ Z, data = grid_df, FUN = function(Z) c(somme = sum(Z), n = length(Z)))
    sum_VegDens <- sum_VegDens %>%
      mutate(PAD = PAD[, "somme"] / (plot_area/(vox_size^2))) # Mean VegDens (total VegDens / voxels count in the area)

    if (debug) {

      ggplot(sum_VegDens, aes(x = Z, y = PAD)) +
        geom_line() +
        geom_point() +
        labs(title = paste0("Sum of ", name, "in function of the height"), x = "Height Z", y = name) +
        theme_minimal()

    }

    profile <- sum_VegDens[, c("Z","PAD")]
    names(profile) <- c("Z", name)

    if (debug) {
      print("profile ok")
    }

    if (export) {
      if (!is.null(output_dir)) {
        write_xlsx(profile, file.path(output_dir, paste0("ForestLayers_1Dprofile_", plot_name,".xlsx")))
      } else {
        message ("Please specify an output directory to export results")
        NULL
      }

    }


  }, error = function(e){
    message("Warning detected during the execution of creates_PAD_profile function", e$message)
    profile <- data.frame()
  })

  return(profile)

}
