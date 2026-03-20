#' Visualization of the full vertical profile
#'
#' @description Creates the visual representation of the vertical layers of your study area.
#'
#' @param finalList R list. Contains the return of final_limits() function.
#' @param col_name character. Name of the value of vegetation density (e.g. PAD, biomass, voxel count). Used for visualization purpose. Default is "VegDens".
#' @param debug logical. For development only.
#' @return Returns a plot with the raw vertical profile of the study area and the fitted Weibull distributions for each layer.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
visu <- function(finalList, col_name = "VegDens", debug = F) {

  tryCatch({

    complet <- bind_rows(finalList, .id = "id") %>% # Stack all the data.frame
      dplyr::mutate(
        id = as.integer(.data$id),
        label = paste("Fitted Weibull - Part", .data$id) # Create legend label
      )

    complet <- complet[order(complet$z), ]

    if (debug) {

    }

    # Layers boundaries (visualization)
    hline_y <- complet %>%
      group_by(.data$id) %>%
      summarise(max_z = max(.data$z, na.rm = TRUE)) %>%  # or min(y)
      pull(max_z)
    hline_y <- hline_y[-length(hline_y)] # Remove the maximum (not a layer limit)

    params_df <- complet %>%
      group_by(.data$id, .data$label) %>%
      summarise(
        z_pos = if (all(is.na(peak_position_weibull))) mean(VegDens_pred) else mean(peak_position_weibull, na.rm = TRUE),
        x_pos = 0.01,
        k = mean(K),
        l = mean(L),
        factor = mean(Scale_Factor),
        offset = mean(Offset),
        .groups = "drop"
      ) %>%
      mutate(param_text = paste0("K = ", round(k,3),
                                 "  L = ", round(l,3),
                                 "  Factor = ", round(factor, 3),
                                 "  Offset = ", round(offset, 3)))

    if (debug) {
      print(hline_y)
      print("offset")
    }

    plot <- ggplot() +
      geom_point(data = complet, aes(x = .data$VegDens, y = z, color = paste0(col_name, " profile")), size = 1) +
      geom_path(data = complet, aes(x = .data$VegDens, y = z, color = paste0(col_name, " profile")), linewidth = 0.5) +
      geom_path(data = complet, aes(x = VegDens_pred, y = z, color = .data$label), linewidth = 0.7) + # Weibull estimations
      geom_hline(aes(yintercept = hline_y, color = "Layer limit"), linetype = "dashed", linewidth = 0.5) +
      geom_point(data = complet, aes(x = VegDens_max_weibull, y = peak_position_weibull, color = "Maximum"), size = 1.5) + # Add maxima
      geom_text(data = params_df, aes(x = .data$x_pos, y = .data$z_pos, label = .data$param_text, color = .data$label),
                hjust = 0, vjust = -0.5, size = 2, show.legend = FALSE) +
      labs(title = paste0("Comparison of ", col_name, " values and fitted Weibull"),
           x = col_name, y = "Height above ground (m)", color = "Legend") +
      theme_minimal() +
      scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
      theme(
        legend.position = c(0.7, 0.5),
        legend.justification = "left",
        plot.margin = margin(5, 15, 5, 5), # top, right, bottom, left
        legend.background = element_rect(fill = "grey95", color = NA),
        legend.key = element_rect(fill = "grey95", color = NA),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.title.x = element_text(size = 8),   # size of the title for the x axe
        axis.title.y = element_text(size = 8),   # size of the title for the Y axe
        axis.text.x  = element_text(size = 5),   # size of the spacing for the X axe
        axis.text.y  = element_text(size = 5),   # size of the spacing for the Y axe
        plot.title   = element_text(size = 9)
      ) +
      guides(colour = guide_legend(ncol = 1)) + # it's possible to display several columns in the legend
      scale_color_manual(
        values = setNames(
          c("gray25", "darkgreen", "chartreuse4", "chartreuse3", "blue", "darkred"),
          c(paste0(col_name, " profile"),
            "Fitted Weibull - Part 1",
            "Fitted Weibull - Part 2",
            "Fitted Weibull - Part 3",
            "Layer limit",
            "Maximum"
          ))
      )

    plot(plot)

  }, error = function(e){
    message("Warning detected during the execution of visu function (creation of the final graph failed)", e$message)
    plot <- list()
  })

  return(plot)

}
