#' Automatic identification of vertical layers in the 1D profile
#'
#' @description Subdivides a VegDens vertical profile into different vegetation layers.
#'
#' @param profile A 2 columns dataframe describing the VegDens vertical distribution (Z, VegDens). Can be the result of the creates_PAD_profile() function.
#' @param diff_threshold numeric. The difference threshold (between 0 and 1) corresponding to a percentage of the segment’s maximum VegDens. Default is 0.25. See future paper for more details.
#' @param diff_max_threshold numeric. The difference threshold (between 0 and 1) corresponding to a percentage of the profile’s maximum VegDens. Default is 0.25. See future paper for more details.
#' @param col_name character. Name of the value of vegetation density (e.g. PAD, biomass, voxel count). Used for visualization purpose. Default is "VegDens".
#' @param print logical. Choose TRUE to visualize intermediate graphs created during the function.
#' @param debug logical. For development only.
#' @param export logical. Choose TRUE and specify a directory (output_dir) to export the 1D profile.
#' @param output_dir character. Output directory to export the intermediate results of the function.
#' @param plot_name character. Name of your data, for proper export procedure.
#' @return Returns a list of subtable of profile. Each subtable is a 2 columns dataframe describing the VegDens vertical distribution for one layer (z, VegDens).
#' @keywords internal
#' @author Amélie Juckler
#'
#'
segments_PAD_profile <- function(profile, diff_threshold = 0.25, diff_max_threshold = 0.25, col_name = "VegDens", print = F, debug = F, export = F, output_dir = NULL, plot_name = NULL){

  tryCatch({

    segmentList = list()

    # Basic verifications

    if (("X" %in% colnames(df)) & ("Y" %in% colnames(df))) {
      print ("Your data is a 3D matrix, use creates_PAD_profile() to convert your matrix into a PAD vertical profile")
    }
    if (!is.data.frame(profile)) stop("Your object must be a 2 columns data.frame.")
    if (ncol(profile) < 2) stop("Column missing.")
    if (any(is.na(profile[, 1:2]))) {
      warning("NA detected and deleted.")
      profile <- na.omit(profile[, 1:2])
    }
    if (debug) {message("profile ok")}

    # Save graphic parameters and function to adjust margins if plot too large

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    adjust_margins <- function(debug = FALSE) {
      try({
        usr_dims <- par("pin")
        if (any(usr_dims < c(3, 3))) {
          par(mar = c(3, 3, 2, 1))
          if (debug) message("Plot margins were too small, adjustment done.")
        } else {
          par(mar = c(4, 4, 2, 1))
        }
      }, silent = TRUE)
    }


    # Step 1 : Identify minima and maxima of the profile

    # Spline
    x = profile[,1]
    y = profile[,2]
    spline_fit <- smooth.spline(x, y, spar = 0.4) # See method in Penner et al. (2023)

    if (debug) {
      adjust_margins(debug)
      try({
        plot(x, y, main = "Smoothing Spline Fit", xlab = "Height (m)", ylab = col_name, pch = 16)
        lines(spline_fit, col = "blue", lwd = 2)
      }, silent = TRUE)
    }

    # First and second derivatives computation
    spline_first_derivative <- predict(spline_fit, x, deriv = 1)
    first_derivative <- spline_first_derivative$y
    spline_second_derivative <- predict(spline_fit, x, deriv = 2)

    # Identification of zeros in the first derivative (min and max)
    zero_crossings <- which(diff(sign(first_derivative)) != 0)

    # Identification of the first min/max
    edge_left_max  <- if (first_derivative[1] < 0) 1 else NULL
    edge_left_min  <- if (first_derivative[1] > 0) 1 else NULL
    filtered_crossings <- zero_crossings

    # Distinction of min and max via second derivative
    # Second derivative < 0 : local maximum ; Seconde derivative > 0 : local minimum
    peaks <- c(edge_left_max, filtered_crossings[spline_second_derivative$y[filtered_crossings] < 0])
    valleys <- c(edge_left_min, filtered_crossings[spline_second_derivative$y[filtered_crossings] > 0])

    # Maxima extraction
    highlight_x_peaks <- spline_first_derivative$x[peaks]
    highlight_y_peaks <- predict(spline_fit, highlight_x_peaks)$y

    # Minima extraction
    highlight_x_valleys <- spline_first_derivative$x[valleys]
    highlight_y_valleys <- predict(spline_fit, highlight_x_valleys)$y

    if (debug) {
      adjust_margins(debug)
      try({
        plot(x, y, main = "Smoothing Spline with min and max",
             xlab = "Height (m)", ylab = col_name, pch = 18)
        lines(spline_fit, col = "blue", lwd = 2)
        points(highlight_x_peaks, highlight_y_peaks, col = "red", pch = 19, cex = 1) # Add maximum in red
        abline(v = highlight_x_peaks, col = "red", lty = 2) # Vertical lines
        points(highlight_x_valleys, highlight_y_valleys, col = "green", pch = 19, cex = 1) # Add minimum in green
        abline(v = highlight_x_valleys, col = "green", lty = 2)
        # Legend
        legend("topright", legend = c("Smoothed spline", "Maxima", "Minima"),
               col = c("blue", "red", "green"), lwd = 2, pch = c(NA, 19, 19), lty = c(1, 2, 2))
      }, silent = TRUE)
    }

    if (debug) {
      print("Affichage des segments initiaux")
    }

    # Step 2 : Remove leftover ground points

    removed_segment <- NULL

    if (length(highlight_x_peaks) > 1 && length(highlight_x_valleys) > 0) {

      first_peak <- min(highlight_x_peaks)
      first_valley <- min(highlight_x_valleys)
      second_valley <- highlight_x_valleys[2]
      diff_height <- abs(first_valley - first_peak)

      if (diff_height < 2 && diff_height > 0) {

        if (debug) message("Leftover ground points detected : height difference between max and min = ", round(diff_height, 2), " m.")

        # Save deleted segment and remove it from the profile

        if (first_peak < first_valley) {
          if (first_valley < 5) { # to avoid the suppression of heigher points

            lower_bound <- first_peak
            upper_bound <- first_valley

            removed_segment <- profile[profile[,1] >= lower_bound & profile[,1] <= upper_bound, ]

            # Delete the segment
            profile <- profile[profile[,1] > upper_bound, ]
            x = profile[,1]
            y = profile[,2]
            highlight_x_peaks <- highlight_x_peaks[-1]
            highlight_y_peaks <- highlight_y_peaks[-1]
            highlight_x_valleys <- highlight_x_valleys[-1]
            highlight_y_valleys <- highlight_y_valleys[-1]


            if (debug) {
              try({
                abline(v = c(first_peak, first_valley), col = "orange", lwd = 2, lty = 3)
                mtext("Deleted segment (<2 m)", side = 3, col = "red")
              }, silent = TRUE)
            }

          }

        } else if (length(second_valley) == 0 || is.na(second_valley)) {

          message("Not enough valleys -> stop")

        } else if (second_valley < 5) {

          lower_bound <- first_valley
          upper_bound <- second_valley

          highlight_x_valleys <- highlight_x_valleys[-1]
          highlight_y_valleys <- highlight_y_valleys[-1]

          removed_segment <- profile[profile[,1] >= lower_bound & profile[,1] <= upper_bound, ]

          # Delete the segment
          profile <- profile[profile[,1] > upper_bound, ]
          x = profile[,1]
          y = profile[,2]
          highlight_x_peaks <- highlight_x_peaks[-1]
          highlight_y_peaks <- highlight_y_peaks[-1]
          highlight_x_valleys <- highlight_x_valleys[-1]
          highlight_y_valleys <- highlight_y_valleys[-1]


          if (debug) {
            try({
              abline(v = c(first_peak, first_valley), col = "orange", lwd = 2, lty = 3)
              mtext("Deleted segment (<2 m)", side = 3, col = "red")
            }, silent = TRUE)
          }

        }

      } else {
        if (debug) message("No anomaly detection (height difference between min and max = ", round(diff_height, 2), " m).")
      }
    } else if (debug) {
      message("Not enough maxima and minima to detect anomalies near the ground.")
    }


    # Step 3 : Compute initial segments

    start_index <- 1

    for (i in seq_along(highlight_x_valleys)) {

      end_index <- which(x == highlight_x_valleys[i])
      segment <- data.frame(x = x[start_index:end_index], y = y[start_index:end_index])

      if (nrow(segment) > 0) {
        segmentList[[length(segmentList) + 1]] <- segment
        start_index <- end_index + 1
        }

    }

    # Add last segment after de last valley
    if (start_index <= length(x)) {
      segment <- data.frame(x = x[start_index:length(x)], y = y[start_index:length(x)])
      if (nrow(segment) > 0) {
        segmentList[[length(segmentList) + 1]] <- segment
      }
    }

    # Compute maxima and difference with the boundaries

    maxima_info <- data.frame(Segment = integer(), Max_X = numeric(), Max_Y = numeric(), Min_X_left = numeric(), Min_Y_left = numeric(), Min_X_right = numeric(), Min_Y_right = numeric(), Diff_left = numeric(), Diff_right = numeric())

    for (i in seq_along(segmentList)) {
      segment <- segmentList[[i]]

      # ignore empty segments
      if (nrow(segment) == 0 || length(segment$x) == 0 || length(segment$y) == 0) {
        next
      }

      max_y <- max(segment$y)
      max_x <- segment$x[which.max(segment$y)]
      min_x_left <- segment$x[1]
      min_x_right <- segment$x[length(segment$x)]
      min_y_left <- segment$y[1]
      min_y_right <- segment$y[length(segment$y)]
      diff_left <- max_y - min_y_left
      diff_right <- max_y - min_y_right
      maxima_info <- rbind(maxima_info, data.frame(Segment = i, Max_X = max_x, Max_Y = max_y, Min_X_left = min_x_left, Min_Y_left = min_y_left, Min_X_right = min_x_right, Min_Y_right = min_y_right, Diff_left = diff_left, Diff_right = diff_right))
    }


    # Visualization of the candidate layers
    if (debug) {
      try({
        if (export) {
          if (!is.null(output_dir)) {

            png(file.path(output_dir, paste0("Initial_segments_", plot_name, ".png")),
                width = 1400, height = 800, res = 150)

          } else {
            message ("Please specify an output directory to export results")
            NULL
          }

        }

        # adjust margins
        par(mar = c(4, 4, 2, 14), xpd = NA)

        plot(
          y, x,
          xlim = c(0, max(y, na.rm = TRUE) + 1),
          ylim = c(0, max(x, na.rm = TRUE) + 1),
          main = "Initial segments",
          xlab = col_name,
          ylab = "Height (m)",
          pch = 16
        )

        sp <- predict(spline_fit, x)
        lines(sp$y, sp$x, col = "blue", lwd = 2)

        nseg <- length(segmentList)
        if (nseg == 0) stop("segmentList is empty, it is not possible to create a plot.")

        colors <- viridis::viridis(nseg)

        for (i in seq_along(segmentList)) {
          seg <- segmentList[[i]]
          if (is.null(seg) || nrow(seg) < 2) {next}
          if (!all(c("x","y") %in% names(seg))) {next}
          if (anyNA(seg$x) || anyNA(seg$y)) {next}
          lines(seg$y, seg$x, col = colors[i], lwd = 2)
        }

        if (!is.null(maxima_info) && nrow(maxima_info) > 0) {
          points(maxima_info$Max_Y, maxima_info$Max_X,
                 col = "red", pch = 19, cex = 1.5)
        }

        # Legend (next to the plot)
        usr <- par("usr")
        x_leg <- usr[2] + 0.03*(usr[2]-usr[1])
        y_leg <- usr[4]
        legend(
          x = x_leg, y = y_leg,
          xjust = 0, yjust = 1,
          legend = c("Smoothed spline", paste("Segment", seq_len(nseg)), "Maxima"),
          col = c("blue", colors, "red"),
          lwd = 2,
          pch = c(NA, rep(NA, nseg), 19),
          lty = c(1, rep(1, nseg), NA),
          bty = "n"
        )

        if (export) dev.off()

      }, silent = TRUE)
    }


    # Step 4 : Reduce to major layers

    segmentList <- segmentList[sapply(segmentList, function(z) nrow(z) > 0)] #remove empty segments

    # 1. Difference with segment boundaries

    repeat {

      nseg <- length(segmentList)

      if (nseg == 0) break

      # Find segments that do not respect the condition
      small_diff <- which(
        (maxima_info$Diff_left  < diff_threshold * maxima_info$Max_Y) |
          (maxima_info$Diff_right < diff_threshold * maxima_info$Max_Y)
      )


      # If the condition is respected -> stop
      if (length(small_diff) == 0) break

      # Start with the first segment
      i <- small_diff[1]

      # Which one is the smallest
      small_left  <- maxima_info$Diff_left[i]  < diff_threshold * maxima_info$Max_Y[i]
      small_right <- maxima_info$Diff_right[i] < diff_threshold * maxima_info$Max_Y[i]

      if (small_left & small_right) {
        if (maxima_info$Diff_left[i] <= maxima_info$Diff_right[i]) {
          small_right <- FALSE
        } else {
          small_left <- FALSE
        }
      }

      # If the smallest is small_left -> combine with the precedent segment
      if (small_left) {
        if (i == 1) {
          # If it's the first segment -> combine with the next segment
          if (nseg >= 2) cible <- 2 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i - 1
        }
        # Limits verification
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid boundaries (left): ", cible)
        }
        # Merge
        # Make sure the heights are sorted in ascending order after the merging step
        segmentList[[cible]] <- rbind(segmentList[[cible]], segmentList[[i]])
      } else if (small_right) {
        if (i == length(segmentList)) {
          # if no right neighbor : merge with left neighbor
          if (nseg >= 2) cible <- i - 1 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i + 1
        }
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid data (right): ", cible)
        }
        segmentList[[cible]] <- rbind(segmentList[[i]], segmentList[[cible]])
      } else {
        # No merging possible (skip)
        warning("No valid neighbor found to merge with the segment", i, " step skipped")
        break
      }

      # Make sure that the heights are well sorted
      segc <- segmentList[[cible]]
      ord <- order(segc$x)
      segmentList[[cible]] <- segc[ord, , drop = FALSE]

      # Delete old segment i (before merging step)
      segmentList <- segmentList[-i]

      # Recomputation of maxima_info

      maxima_info <- data.frame(
        Segment = integer(),
        Max_X = numeric(),
        Max_Y = numeric(),
        Min_X_left = numeric(),
        Min_Y_left = numeric(),
        Min_X_right = numeric(),
        Min_Y_right = numeric(),
        Diff_left = numeric(),
        Diff_right = numeric()
      )

      for (j in seq_along(segmentList)) {
        seg <- segmentList[[j]]

        max_y <- max(seg$y)
        max_x <- seg$x[which.max(seg$y)]

        min_x_left  <- seg$x[1]
        min_y_left  <- seg$y[1]
        min_x_right <- seg$x[length(seg$x)]
        min_y_right <- seg$y[length(seg$y)]

        diff_left  <- max_y - min_y_left
        diff_right <- max_y - min_y_right

        maxima_info <- rbind(
          maxima_info,
          data.frame(
            Segment = j,
            Max_X = max_x,
            Max_Y = max_y,
            Min_X_left = min_x_left,
            Min_Y_left = min_y_left,
            Min_X_right = min_x_right,
            Min_Y_right = min_y_right,
            Diff_left = diff_left,
            Diff_right = diff_right
          )
        )
      }
    } # repeat loop 1


    # 2. Vertical extent

    repeat {

      nseg <- length(segmentList)

      if (nseg == 0) break

      # Find segments that do not respect the condition
      small_height_range <- which(abs(maxima_info$Min_X_right - maxima_info$Min_X_left) < 3)

      # If the condition is respected -> stop
      if (length(small_height_range) == 0) break

      # Start with the first segment
      i <- small_height_range[1]

      # Which one is the smallest
      small_left  <- TRUE
      small_right <- TRUE

      if (small_left & small_right) {
        if (maxima_info$Diff_left[i] <= maxima_info$Diff_right[i]) {
          small_right <- FALSE
        } else {
          small_left <- FALSE
        }
      }

      # If the smallest is small_left -> combine with the precedent segment
      if (small_left) {
        if (i == 1) {
          # If it's the first segment -> combine with the next segment
          if (nseg >= 2) cible <- 2 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i - 1
        }
        # Limits verification
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid boundaries (left): ", cible)
        }
        # Merge
        # Make sure the heights are sorted in ascending order after the merging step
        segmentList[[cible]] <- rbind(segmentList[[cible]], segmentList[[i]])
      } else if (small_right) {
        if (i == length(segmentList)) {
          # if no right neighbor : merge with left neighbor
          if (nseg >= 2) cible <- i - 1 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i + 1
        }
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid data (right): ", cible)
        }
        segmentList[[cible]] <- rbind(segmentList[[i]], segmentList[[cible]])
      } else {
        # No merging possible (skip)
        warning("No valid neighbor found to merge with the segment ", i, "; step skipped")
        break
      }

      # Make sure that the heights are well sorted
      segc <- segmentList[[cible]]
      ord <- order(segc$x)
      segmentList[[cible]] <- segc[ord, , drop = FALSE]

      # Delete old segment i
      segmentList <- segmentList[-i]

      # Recomputation of maxima_info

      maxima_info <- data.frame(
        Segment = integer(),
        Max_X = numeric(),
        Max_Y = numeric(),
        Min_X_left = numeric(),
        Min_Y_left = numeric(),
        Min_X_right = numeric(),
        Min_Y_right = numeric(),
        Diff_left = numeric(),
        Diff_right = numeric()
      )

      for (j in seq_along(segmentList)) {
        seg <- segmentList[[j]]

        max_y <- max(seg$y)
        max_x <- seg$x[which.max(seg$y)]

        min_x_left  <- seg$x[1]
        min_y_left  <- seg$y[1]
        min_x_right <- seg$x[length(seg$x)]
        min_y_right <- seg$y[length(seg$y)]

        diff_left  <- max_y - min_y_left
        diff_right <- max_y - min_y_right

        maxima_info <- rbind(
          maxima_info,
          data.frame(
            Segment = j,
            Max_X = max_x,
            Max_Y = max_y,
            Min_X_left = min_x_left,
            Min_Y_left = min_y_left,
            Min_X_right = min_x_right,
            Min_Y_right = min_y_right,
            Diff_left = diff_left,
            Diff_right = diff_right
          )
        )
      }
    } # repeat loop 2

    # 3. Difference between max of the segment and max of the profile

    max_profile <- max(maxima_info$Max_Y)

    repeat {

      nseg <- length(segmentList)

      if (debug) {
        print(maxima_info$Max_Y)
        print(max_profile)
      }

      if (nseg == 0) break

      # Find segments that do not respect the condition
      small_diff_y <- which(maxima_info$Max_Y < (diff_max_threshold * max_profile))

      # If the condition is respected -> stop
      if (length(small_diff_y) == 0) break

      # Start with the first segment
      i <- small_diff_y[1]

      # Which one is the smallest
      small_left  <- TRUE
      small_right <- TRUE

      if (small_left & small_right) {
        if (maxima_info$Diff_left[i] <= maxima_info$Diff_right[i]) {
          small_right <- FALSE
        } else {
          small_left <- FALSE
        }
      }

      # If the smallest is small_left -> combine with the precedent segment
      if (small_left) {
        if (i == 1) {
          # If it's the first segment -> combine with the next segment
          if (nseg >= 2) cible <- 2 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i - 1
        }
        # Limits verification
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Problem with the neighbor segment (left): ", cible)
        }
        # Merge
        # Make sure that the heights are sorted in ascending order after the merging step
        segmentList[[cible]] <- rbind(segmentList[[cible]], segmentList[[i]])
      } else if (small_right) {
        if (i == length(segmentList)) {
          # if no right neighbor : merge with left neighbor
          if (nseg >= 2) cible <- i - 1 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i + 1
        }
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid data (right): ", cible)
        }
        segmentList[[cible]] <- rbind(segmentList[[i]], segmentList[[cible]])
      } else {
        # No merging possible (skip)
        warning("No valid neighbor found to merge with the segment ", i, "; step skipped")
        break
      }

      # Make sure that the heights are well sorted
      segc <- segmentList[[cible]]
      ord <- order(segc$x)
      segmentList[[cible]] <- segc[ord, , drop = FALSE]

      # Delete old segment i
      segmentList <- segmentList[-i]

      # Recomputation of maxima_info

      maxima_info <- data.frame(
        Segment = integer(),
        Max_X = numeric(),
        Max_Y = numeric(),
        Min_X_left = numeric(),
        Min_Y_left = numeric(),
        Min_X_right = numeric(),
        Min_Y_right = numeric(),
        Diff_left = numeric(),
        Diff_right = numeric()
      )

      for (j in seq_along(segmentList)) {
        seg <- segmentList[[j]]

        max_y <- max(seg$y)
        max_x <- seg$x[which.max(seg$y)]

        min_x_left  <- seg$x[1]
        min_y_left  <- seg$y[1]
        min_x_right <- seg$x[length(seg$x)]
        min_y_right <- seg$y[length(seg$y)]

        diff_left  <- max_y - min_y_left
        diff_right <- max_y - min_y_right

        maxima_info <- rbind(
          maxima_info,
          data.frame(
            Segment = j,
            Max_X = max_x,
            Max_Y = max_y,
            Min_X_left = min_x_left,
            Min_Y_left = min_y_left,
            Min_X_right = min_x_right,
            Min_Y_right = min_y_right,
            Diff_left = diff_left,
            Diff_right = diff_right
          )
        )
      }
    } # repeat loop 3

    # 4. Combine segments for which the position (Z) of the peaks are too close

    repeat {

      nseg <- length(segmentList)

      if (nseg == 0) break

      # Check the condition

      # Initialization
      peak_issue_left <- rep(FALSE, nrow(maxima_info))
      peak_issue_right <- rep(FALSE, nrow(maxima_info))

      for (i in seq_len(nrow(maxima_info))) {

        peak_i <- maxima_info$Max_X[i]

        # left neighbor
        if (i > 1) {
          peak_left <- maxima_info$Max_X[i - 1]
          if (abs(peak_i - peak_left) < 5) {
            peak_issue_left[i] <- TRUE
          }
        }

        # right neighbor
        if (i < nrow(maxima_info)) {
          peak_right <- maxima_info$Max_X[i + 1]
          if (abs(peak_right - peak_i) < 5) {
            peak_issue_right[i] <- TRUE
          }
        }
      }

      # Find segments that do not respect the condition
      close_to_peak <- which(peak_issue_left | peak_issue_right)

      # If the condition is respected -> stop
      if (length(close_to_peak) == 0) break

      # Start with the first segment
      i <- close_to_peak[1]

      # Which one is problematic
      fuse_left  <- peak_issue_left[i]
      fuse_right <- peak_issue_right[i]

      if (fuse_left & fuse_right) {

        dist_left  <- abs(maxima_info$Max_X[i] - maxima_info$Max_X[i - 1])
        dist_right <- abs(maxima_info$Max_X[i] - maxima_info$Max_X[i + 1])

        if (dist_left <= dist_right) {
          fuse_right <- FALSE
        } else {
          fuse_left <- FALSE
        }
      }

      # If the smallest is fuse_left -> combine with the precedent segment
      if (fuse_left) {

        if (i == 1) {
          # If it's the first segment -> combine with the next segment
          if (nseg >= 2) cible <- 2 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i - 1
        }

        # Limits verification
        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid boundaries (left): ", cible)
        }

        # Merge
        segmentList[[cible]] <- rbind(segmentList[[cible]], segmentList[[i]])

      } else if (fuse_right) {

        if (i == length(segmentList)) {
          # no right neigbor : fuse with the left one
          if (nseg >= 2) cible <- i - 1 else { warning("Impossible to merge a unique segment"); break }
        } else {
          cible <- i + 1
        }

        if (!(cible %in% seq_len(length(segmentList)))) {
          stop("Invalid data (right): ", cible)
        }

        segmentList[[cible]] <- rbind(segmentList[[i]], segmentList[[cible]])

      } else {

        # No merging possible (skip)
        warning("No valid neighbor found to merge with the segment ", i, "; step skipped")
        break

      }

      # Make sure that the heights are well sorted
      segc <- segmentList[[cible]]
      ord <- order(segc$x)
      segmentList[[cible]] <- segc[ord, , drop = FALSE]

      # Delete old segment i
      segmentList <- segmentList[-i]

      # Recomputation of maxima_info

      maxima_info <- data.frame(
        Segment = integer(),
        Max_X = numeric(),
        Max_Y = numeric(),
        Min_X_left = numeric(),
        Min_Y_left = numeric(),
        Min_X_right = numeric(),
        Min_Y_right = numeric(),
        Diff_left = numeric(),
        Diff_right = numeric()
      )

      for (j in seq_along(segmentList)) {
        seg <- segmentList[[j]]

        max_y <- max(seg$y)
        max_x <- seg$x[which.max(seg$y)]

        min_x_left  <- seg$x[1]
        min_y_left  <- seg$y[1]
        min_x_right <- seg$x[length(seg$x)]
        min_y_right <- seg$y[length(seg$y)]

        diff_left  <- max_y - min_y_left
        diff_right <- max_y - min_y_right

        maxima_info <- rbind(
          maxima_info,
          data.frame(
            Segment = j,
            Max_X = max_x,
            Max_Y = max_y,
            Min_X_left = min_x_left,
            Min_Y_left = min_y_left,
            Min_X_right = min_x_right,
            Min_Y_right = min_y_right,
            Diff_left = diff_left,
            Diff_right = diff_right
          )
        )
      }
    } # repeat loop 4

    # Visualization of the effective layers
    if (print) {

      try({
        print(maxima_info)

        if (export) {
          if (!is.null(output_dir)) {

            png(
              file.path(output_dir, paste0("Intermediate_segments_", plot_name, ".png")),
              width = 1400, height = 800, res = 150)

          } else {
            message ("Please specify an output directory to export results")
            NULL
          }

        }

        # Reverse the axes : X = VegDens, y = height (z)
        par(mar = c(4, 4, 2, 14), xpd = NA)

        plot(
          y, x,
          main = "Effective layers",
          xlab = col_name,
          ylab = "Height (m)",
          pch = 16
        )

        sp <- predict(spline_fit, x)
        lines(sp$y, sp$x, col = "blue", lwd = 2)

        # Define colors scale
        nseg <- length(segmentList)
        colors <- viridis::viridis(max(1, nseg))

        for (i in seq_along(segmentList)) {
          seg <- segmentList[[i]]
          if (is.null(seg) || nrow(seg) < 2) {next}
          if (!all(c("x","y") %in% names(seg))) {next}
          if (anyNA(seg$x) || anyNA(seg$y)) {next}
          lines(seg$y, seg$x, col = colors[i], lwd = 2)
        }

        if (!is.null(maxima_info) && nrow(maxima_info) > 0) {
          points(maxima_info$Max_Y, maxima_info$Max_X, col = "red", pch = 19)
        }

        # Legend next to the plot (avoid overlapping)
        usr <- par("usr")
        x_leg <- usr[2] + 0.03 * (usr[2] - usr[1])
        y_leg <- usr[4]

        legend(
          x = x_leg, y = y_leg,
          xjust = 0, yjust = 1,
          legend = c("Smoothed spline", paste("Segment", seq_len(nseg)), "Maxima"),
          col = c("blue", colors[seq_len(nseg)], "red"),
          lwd = 2,
          pch = c(NA, rep(NA, nseg), 19),
          lty = c(1, rep(1, nseg), NA),
          bty = "n"
        )

        if (export) {
          dev.off()
        }

      }, silent = TRUE)
    }

    gc()

  }, error = function(e){
    message("Warning detected during the execution of segment_PAD_profile function", e$message)
    segmentList <- list()
  })

  return(list(
    segmentList,
    removed_segment))

}
