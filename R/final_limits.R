#' Get final boundaries of the vertical layers
#'
#' @description Redefines each vegetation layer to match the Weibull estimations.
#' weibullList is a list of 25 columns data.frame.
#'
#' @param weibullList R list. Contains the return of weibull_fit() function.
#' @param col_name character. Name of the value of vegetation density (e.g. PAD, biomass, voxel count). Used for visualization purpose. Default is "VegDens".
#' @param print logical. Choose TRUE to visualize intermediate graphs created during the function.
#' @param debug logical. For development only.
#' @param export logical. Choose TRUE and specify a directory (output_dir) to export the 1D profile.
#' @param output_dir character. Output directory to export the intermediate results of the function.
#' @param plot_name character. Name of your data, for proper export procedure.
#' @return Returns a list of 25 columns data.frames with the VegDens estimations, structure parameters and final layers limits.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
final_limits <- function(weibullList, col_name = "VegDens", print = F, debug = F, export = F, output_dir = NULL, plot_name = NULL) {

  tryCatch({

    finalList <- list()

    if (length(weibullList) > 1) {

      # Find Weibull intersection
      find_intersection <- function(param1, param2, x_range) {
        df <- data.frame(x = seq(x_range[1], x_range[2], length.out = 1000))
        df$y1 <- my_weibull(param1, df$x)
        df$y2 <- my_weibull(param2, df$x)

        idx <- which(diff(sign(df$y1 - df$y2)) != 0)

        if (length(idx) == 0) return(NULL)  # There is no crossing

        x_cross <- df$x[idx[1]]
        return(x_cross)
      }

      # Extend weibull estimations to the crossing point
      get_params <- function(params) {

        K <- params$K[1]
        L  <- params$L[1]
        scale_factor  <- params$Scale_Factor[1]
        offset  <- params$Offset[1]

        parametres <- c(K,L,scale_factor,offset)
        return(parametres)

      }

      # Check if crossing is in the delimited range and adjust layers limits
      for (i in 1:(length(weibullList) - 1)) {
        seg1 <- as.data.frame(weibullList[i])
        seg2 <- as.data.frame(weibullList[i + 1])

        # Extract Weibull parameters
        param1 <- get_params(seg1)
        param2 <- get_params(seg2)

        # Extract the boundaries x of the segments
        x_limit <- max(seg1$z, na.rm = TRUE)
        x_max2 <- max(seg2$z, na.rm = TRUE)
        x_min1 <- min(seg1$z, na.rm = TRUE)

        # Security check
        if (!is.finite(x_limit)) {
          warning(paste("No initial limit found between the segments", seg1, "and", seg2))
          final_limits <- c(final_limits, NA)
          next
        }

        if (debug) {
          cat("Initial limit : ", x_limit, "\n")
        }

        # Computing the searching limit for the crossing point
        segment_size1 <- x_limit - x_min1
        segment_size2 <- x_max2 - x_limit

        threshold1 <- segment_size1 * 0.3 # 30% of the vertical extent of the segment
        threshold2 <- segment_size2 * 0.3

        limit1 <- x_limit - threshold1
        limit2 <- x_limit + threshold2

        x_range <- c(limit1, limit2) # Searching limits

        # Updating segments boundaries and VegDens estimation (Weibull) to cover the searching zone
        seg <- rbind(seg1, seg2)
        x1 <- data.frame(x = seg$z[seg$z >= x_min1 & seg$z <= limit2])
        x2 <- data.frame(x = seg$z[seg$z >= limit1 & seg$z <= x_max2])

        if (debug) {
          print(paste("limit1 = ", limit1))
          print(paste("limit2 = ", limit2))
          print(x1)
        }

        curv1 <- my_weibull(param1, x1)
        curve1 <- cbind(x1, curv1)
        colnames(curve1) <- c("z", "VegDens_pred")

        curv2 <- my_weibull(param2, x2)
        curve2 <- cbind(x2, curv2)
        colnames(curve2) <- c("z", "VegDens_pred")


        # Finding crossing point (if present in the searching range)
        x_cross <- find_intersection(param1, param2, x_range)

        if (debug) {
          cat("segment_size1", segment_size1, "\n")
          cat("x_range : ", x_range, "\n")
          cat("x_cross : ", x_cross, "\n")
        }

        # Visualization
        if (print) {

          seg <- seg[order(seg$z), ]

          if (!is.null(x_cross)) {

            plot <- ggplot() +
              geom_point(data = seg, aes(x = .data$VegDens, y = z, color = col_name), size = 1) +
              geom_path(data = seg, aes(x = .data$VegDens, y = z, color = col_name), linewidth = 0.5) +
              geom_path(data = curve2, aes(x = .data$VegDens_pred, y = z, color = "Fitted Weibull - part 2"), linewidth = 0.7) +
              geom_path(data = curve1, aes(x = .data$VegDens_pred, y = z, color = "Fitted Weibull - part 1"), linewidth = 0.7) +  # Courbe prédite
              geom_hline(aes(yintercept = c(x_range[1], x_range[2]), color = "Searching limits"), linetype = "dashed", linewidth = 0.5) +
              geom_hline(aes(yintercept = x_cross, color = "Crossing point"), linetype = "dashed", linewidth = 0.5) +
              labs(title = paste0("Comparison of ", col_name, " values and fitted Weibull"),
                   x = col_name, y = "Height above ground (m)", color = "Legend") +
              theme_minimal() +
              theme(legend.position = c(0.75, 0.25)) +
              scale_color_manual(
                values = c(
                  col_name = "gray25",
                  "Fitted Weibull - part 1" = "darkgreen",
                  "Fitted Weibull - part 2" = "chartreuse4",
                  "Searching limits" = "dodgerblue3",
                  "Crossing point" = "blue"
                )
              )

            plot(plot)

            if (export) {
              if (!is.null(output_dir)) {
                ggsave(file.path(output_dir, paste0("Final_limits_", plot_name, "_part_", i, ".png")), plot)
              } else {
                message ("Please specify an output directory to export results")
                NULL
              }

            }

          } else {

            plot <- ggplot() +
              geom_point(data = seg, aes(x = .data$VegDens, y = z, color = col_name)) +
              geom_path(data = seg, aes(x = .data$VegDens, y = z, color = col_name), linewidth = 0.5) +
              geom_path(data = curve2, aes(x = .data$VegDens_pred, y = z, color = "Fitted Weibull - part 2"), linewidth = 0.7) +
              geom_path(data = curve1, aes(x = .data$VegDens_pred, y = z, color = "Fitted Weibull - part 1"), linewidth = 0.7) +  # Courbe prédite
              geom_hline(aes(yintercept = c(x_range[1], x_range[2]), color = "Searching limits"), linetype = "dashed", linewidth = 0.5) +
              geom_hline(aes(yintercept = x_limit, color = "Initial limit (no update)"), linetype = "dashed", linewidth = 0.5) +
              labs(title = paste0("Comparison of ", col_name, " values and fitted Weibull"),
                   x = col_name, y = "Height above ground (m)", color = "Legend") +
              theme_minimal() +
              theme(legend.position = c(0.75, 0.25)) +  # Coordonnées (x, y) de la légende dans la zone de tracé
              scale_color_manual(
                values = c(
                  col_name = "gray25",
                  "Fitted Weibull - part 1" = "darkgreen",
                  "Fitted Weibull - part 2" = "chartreuse4",
                  "Searching limits" = "dodgerblue3",
                  "Initial limit (no update)" = "blue"
                )
              )

            plot(plot)

            if (export) {
              if (!is.null(output_dir)) {
                ggsave(file.path(output_dir, paste0("Final_limits_", plot_name, "_part_", i, ".png")), plot)
              } else {
                message ("Please specify an output directory to export results")
                NULL
              }
            }

          }

        }

        # Create the final vector containing every layer boundaries
        if ((!is.null(x_cross)) & (seg1$NRMSE[1] < 0.15) & (seg2$NRMSE[1] < 0.15)) {

          if (i == 1) {
            final_limits <- c(min(seg1$z, na.rm = TRUE), x_cross)
          } else {
            final_limits <- c(final_limits, x_cross)
          }

          if (i == length(weibullList) - 1) {
            final_limits <- c(final_limits, max(seg2$z, na.rm = TRUE))
          }

        } else {

          # Keep the initial limit
          if (i == 1) {
            final_limits <- c(min(seg1$z, na.rm = TRUE), x_limit)
          } else {
            final_limits <- c(final_limits, x_limit)
          }

          if (i == length(weibullList) - 1) {
            final_limits <- c(final_limits, max(seg2$z, na.rm = TRUE))
          }

        }

        if (debug) {
          print(paste("Final vector, layers boundaries : ", final_limits))
        }

      } # end of for loop

      for (i in 1:(length(weibullList))) {

        seg <- do.call(rbind, weibullList)
        seg1 <- data.frame(weibullList[i])
        x1 <- seg[seg$z >= final_limits[i] & seg$z <= final_limits[i+1],]

        if (debug) {
          print(paste("limit1 = ", final_limits[i]))
          print(paste("limit2 = ", final_limits[i+1]))
        }

        param1 <- param1 <- get_params(seg1)
        curv1 <- data.frame(my_weibull(param1, x1$z))
        curve1 <- cbind(x1$z, x1$VegDens, curv1)

        colnames(curve1) <- c("z", "VegDens", "VegDens_pred")


        # Metrics computation
        calculate_area <- function(x, y) {
          dx <- diff(x)
          avg_height <- (y[-1] + y[-length(y)]) / 2
          abs(sum(dx * avg_height))
        }

        nrmse <- (sqrt(mean((curve1$VegDens - curve1$VegDens_pred)^2, na.rm = TRUE)))/(max(curve1$VegDens) - min(curve1$VegDens))
        VegDens_max_init <- max(curve1$VegDens)
        peak_position_init <- curve1$z[which.max(curve1$VegDens)]
        total_VegDens_init <- calculate_area(curve1$z, curve1$VegDens)
        total_VegDens_weibull <- calculate_area(curve1$z, curve1$VegDens_pred)
        VegDens_mean_init <- mean(curve1$VegDens, na.rm = TRUE)
        VegDens_mean_weibull <- mean(curve1$VegDens_pred, na.rm = TRUE)
        VegDens_med_init <- median(curve1$VegDens, na.rm = TRUE)
        VegDens_med_weibull <- median(curve1$VegDens_pred, na.rm = TRUE)
        VegDens_var_init <- var(curve1$VegDens, na.rm = TRUE)
        VegDens_var_weibull <- var(curve1$VegDens_pred, na.rm = TRUE)
        quantiles <- quantile(curve1$VegDens, probs = c(0.05,0.25,0.75,0.95), na.rm = TRUE, names = FALSE, type = 7)
        skewness <- (VegDens_mean_init - VegDens_med_init) / sqrt(VegDens_var_init) # Pearson coefficient (median)

        if (debug) {
          print(paste("NRMSE : ", nrmse))
        }

        segment_final <- data.frame(
          z = curve1$z,
          VegDens = curve1$VegDens,
          VegDens_pred = curve1$VegDens_pred,
          K = param1[1],
          L = param1[2],
          Scale_Factor = param1[3],
          Offset = param1[4],
          VegDens_max_init = VegDens_max_init,
          VegDens_max_weibull = max(curve1$VegDens_pred),
          peak_position_init = peak_position_init,
          peak_position_weibull = curve1$z[which.max(curve1$VegDens_pred)],
          NRMSE = nrmse,
          total_VegDens_init = total_VegDens_init,
          total_VegDens_weibull = total_VegDens_weibull,
          VegDens_mean_init = VegDens_mean_init,
          VegDens_mean_weibull = VegDens_mean_weibull,
          VegDens_med_init = VegDens_med_init,
          VegDens_med_weibull = VegDens_med_weibull,
          VegDens_var_init = VegDens_var_init,
          VegDens_var_weibull = VegDens_var_weibull,
          quantile_05 = quantiles[1],
          quantile_25 = quantiles[2],
          quantile_75 = quantiles[3],
          quantile_95 = quantiles[4],
          skewness = skewness
        )

        finalList[[i]] <- segment_final

      }

    } else {

      finalList <- weibullList

    }

    if (debug) {

    }

  }, error = function(e){
    message("Warning detected during the execution of final_limits function", e$message)
    finalList <- list()
  })

  return(finalList)

}
