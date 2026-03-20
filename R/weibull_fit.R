#' Fits a Weibull distribution on each vertical layer
#'
#' @description Fits Weibull distribution on the VegDens profile in each vegetation layer.
#' segmentList is a list of 2 columns data.frame corresponding to the VegDens profile of different vegetation layers (x = Z, y = VegDens).
#'
#' @param segmentList R list. Contains the return of segments_PAD_profile() function.
#' @param col_name character. Name of the value of vegetation density (e.g. PAD, biomass, voxel count). Used for visualization purpose. Default is "VegDens".
#' @param print logical. Choose TRUE to visualize intermediate graphs created during the function.
#' @param debug logical. For development only.
#' @return Returns a list of subtable of profile. Each subtable is a 25 columns data.frame describing the VegDens vertical distribution, the Weibull distribution and the structural parameters for one layer.
#' @keywords internal
#' @author Amélie Juckler
#'
#'
weibull_fit <- function(segmentList, col_name = "VegDens", print = FALSE, debug = FALSE) {

  tryCatch({

    weibullList <- list()

    plots <- list()

    removed_segment <- segmentList[[2]]
    segmentListb <- segmentList[[1]]
    i <- 1

    for (i in seq_along(segmentListb)) {

      fitted_results <- data.frame()

      data <- segmentListb[[i]]

      x <- data[[1]]
      y <- data[[2]]
      max_x <- data$x[which.max(data$y)]

      if (debug) {
        print("Maximum is ")
        print(max_x)
      }


      # Initialization of the Weibull parameters

      if (length(segmentListb) == 1) {

        init_kl_array <- c(min(x) + 3, max_x - 0.1, max(y) - min(y), min(y) - 1) # other option : c(min(x)+2, data$x[which.max(data$y)], max(y) - min(y), min(y)-1)

      } else {

        init_kl_array <- c(min(x) + 5, max_x - 0.1, max(y) - min(y), min(y) - 1)

      }

      nrmse <- 1
      attempt <- 1

      while (attempt <= 7 && nrmse > 0.2) {

        if (attempt == 1) {
          par_init <- init_kl_array
        } else {
          par_init <- init_kl_array
          par_init[1] <- init_kl_array[1] - 1
          par_init[2] <- init_kl_array[2] - 0.001
        }

        if (debug) {
          print(paste("Peak position (x), segment ", i))
          print(max_x)
        }

        # Find the best Weibull parameters
        residual <- function(params, x, y) {
          predicted <- my_weibull(params, x)
          return(y - predicted)
        }

        fit_result <- tryCatch({
          nls.lm(
            par = init_kl_array, fn = residual, x = x, y = y, lower = c(1.01, 1.01, min(y), min(y)), upper = c(Inf, Inf, Inf, Inf))
        }, error = function(e) {
          if (debug) message("Error, attempt number ", attempt, " : ", e$message)
          NULL
        })

        # If successfull Weibull fitting :
        if (!is.null(fit_result)) {
          fitted_kl_array <- fit_result$par

          if (debug) {
            print(paste("Parameters of the segment ", i, " : "))
            print(fitted_kl_array)
          }

          # Parameters saving
          fitted_k <- fitted_kl_array[1]
          fitted_l <- fitted_kl_array[2]
          fitted_scale_factor <- fitted_kl_array[3]
          fitted_offset <- fitted_kl_array[4]

          # Computation of estimated VegDens (Weibull)
          # Start at x = 0 if first layer
          if (i == 1) {

            predicted_x <- seq(0, round(max(x),1), by = 0.1)
            predicted_y <- my_weibull(fitted_kl_array, predicted_x)
            xf <- c(removed_segment[,1],x)
            yf <- c(removed_segment[,2],y)

          } else {

            predicted_x <- seq(round(min(x),1), round(max(x),1), by = 0.1)
            predicted_y <- my_weibull(fitted_kl_array, predicted_x)
            xf <- x
            yf <- y

          }

          # Interpolation of the estimations for actual x (z)
          y_pred_interp <- approx(predicted_x, predicted_y, xout = xf, rule = 2)$y # rule approximates the values outside the boundaries

          # NRMSE
          nrmse <- (sqrt(mean((yf - y_pred_interp)^2, na.rm = TRUE)))/(max(yf) - min(yf))

          if (attempt == 1) {
            first_fit <- fitted_kl_array
            first_nrmse <- nrmse
            first_predicted_x <- predicted_x
            first_predicted_y <- predicted_y
            first_y_pred_interp <- y_pred_interp
          }

          if (debug) cat("NRMSE = ", round(nrmse, 4), "\n")
        } else {
          if (debug) message("Adjustment failed on attempt ", attempt)
        }

        # If NRMSE is too high (> 0.2), try again
        if (nrmse > 0.2) {
          attempt <- attempt + 1
        } else {
          break
        }
      }

      # If NRMSE id still too high, restores the parameters of the first attempt
      if (nrmse > 0.2) {
        message("NRMSE is still > 0.2 after ", attempt - 1, " attempts. Restores the parameters of the first attempt.")

        nrmse <- first_nrmse
        fitted_kl_array <- first_fit
        predicted_x <- first_predicted_x
        predicted_y <- first_predicted_y
        y_pred_interp <- first_y_pred_interp
      }

      # Metrics computation
      calculate_area <- function(x, y) {
        dx <- diff(x)
        avg_height <- (y[-1] + y[-length(y)]) / 2
        abs(sum(dx * avg_height))
      }

      fitted_results <- data.frame(
        z = xf,
        VegDens = yf,
        VegDens_pred = y_pred_interp,
        K = fitted_kl_array[1],
        L = fitted_kl_array[2],
        Scale_Factor = fitted_kl_array[3],
        Offset = fitted_kl_array[4],
        VegDens_max_init = max(y),
        VegDens_max_weibull = max(y_pred_interp),
        peak_position_init = x[which.max(y)],
        peak_position_weibull = xf[which.max(y_pred_interp)],
        NRMSE = nrmse,
        total_VegDens_init <- calculate_area(x, y),
        total_VegDens_weibull <- calculate_area(x, y_pred_interp),
        VegDens_mean_init <- mean(y, na.rm = TRUE),
        VegDens_mean_weibull <- mean(y_pred_interp, na.rm = TRUE),
        VegDens_med_init <- median(y, na.rm = TRUE),
        VegDens_med_weibull <- median(y_pred_interp, na.rm = TRUE),
        VegDens_var_init <- var(y, na.rm = TRUE),
        VegDens_var_weibull <- var(y_pred_interp, na.rm = TRUE),
        quantile_05 <- quantile(y, probs = c(0.05), na.rm = TRUE, names = FALSE, type = 7),
        quantile_25 <- quantile(y, probs = c(0.25), na.rm = TRUE, names = FALSE, type = 7),
        quantile_75 <- quantile(y, probs = c(0.75), na.rm = TRUE, names = FALSE, type = 7),
        quantile_95 <- quantile(y, probs = c(0.95), na.rm = TRUE, names = FALSE, type = 7),
        skewness <- (VegDens_mean_init - VegDens_med_init) / sqrt(VegDens_var_init) # Pearson coefficient (median)
      )


      weibullList[[i]] <- fitted_results


      # Visualization
      if (print) {

        plot <- ggplot() +
          geom_point(aes(x = xf, y = yf), color = 'blue') +
          geom_line(aes(x = xf, y = yf), color = 'blue') +
          geom_line(aes(x = predicted_x, y = predicted_y), color = 'green') +
          geom_point(aes(x = max_x, y = max(data$y)), color = 'red') +
          labs(title = paste("Weibull fitting - Part", i),
               x = "z", y = col_name) +
          theme_minimal()

        print(plot)

      }
    }

    # Global NRMSE
    complet <- do.call(rbind, weibullList)

    if (debug) {

    }

    global_nrmse <- (sqrt(mean((complet$VegDens - complet$VegDens_pred)^2, na.rm = TRUE)))/(max(complet$VegDens) - min(complet$VegDens)) #PAD - PAD_pred
    print(paste("Global NRMSE : ", global_nrmse))

    if (print) {

      complet_ <- complet[order(complet$z), ]

      weib <- ggplot(complet_, aes(x = .data$VegDens, y = z)) +
        geom_point(aes(color = col_name)) +
        geom_path(aes(color = col_name), group = 1) +
        geom_path(aes(x = VegDens_pred, y = z, color = "Fitted Weibull"), group = 1) +  # Estimated VegDens
        labs(title = paste0("Comparison of ", col_name, "values and fitted Weibull"),
             x = col_name, y = "Height above ground (m)", color = "Legend") +
        theme_minimal() +
        theme(legend.position = c(0.75, 0.25)) + # Legend box coordinates
        geom_point(aes(x = VegDens_max_weibull, y = peak_position_weibull, color = "Maxima"), size = 2)
      weib <- weib +
        scale_color_manual(
          values = setNames(
            c("darkgreen", "green", "red"),
            c(col_name, "Fitted Weibull", "Maxima")
          )
        )

      plot(weib)
    }

  }, error = function(e){
    message("Warning detected during the execution of weibull_fit function", e$message)
    weibullList <- list()
    nrmse <- NA
  })

  return(list(
    weibullList = weibullList,
    nrmse = global_nrmse
  ))

}
