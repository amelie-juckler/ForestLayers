test_that("VerticalLayers_from_1Dprofile works well", {

  profile1 <- data.frame(
    X = 1:100,
    VegDens = runif(100, 0, 1)
  )

  # Define ForestLayers parameters
  diff_threshold = 0.25
  diff_max_threshold = 0.25
  col_name <- "VegDens (units)"
  output_dir <- tempdir()
  plot_name <- "my_study_area"

  # Function

  png(tempfile(fileext = ".png"))
  on.exit(dev.off(), add = TRUE)

  my_layers <- VerticalLayers_from_1Dprofile(
    profile = profile1,
    diff_threshold,
    diff_max_threshold,
    col_name,
    print = FALSE,
    debug = FALSE,
    export = TRUE,
    output_dir,
    plot_name
  )

  # Check output
  expect_type(my_layers, "list")
  expect_named(my_layers, c("finalResults", "nrmse", "r2", "final_plot"))
  expect_true(is.numeric(my_layers$nrmse))
  expect_true(is.numeric(my_layers$r2))
  expect_s3_class(my_layers$final_plot, "ggplot")

  # Export function
  results <- export_DVV(ForestLayersList = my_layers, output_dir, plot_name)

  # Check output
  expect_true(file.exists(file.path(output_dir, paste0("ForestLayers_", plot_name, ".xlsx"))))
  expect_gt(file.info(file.path(output_dir, paste0("ForestLayers_", plot_name, ".xlsx")))$size, 0)

  # Check png file
  png_path <- file.path(output_dir, "ForestLayers_visualization_my_study_area.png")

  expect_true(file.exists(png_path))

  if (file.info(png_path)$size <= 0) {
    stop("Empty PNG file")
  }
  expect_true(file.info(png_path)$size > 0)

})

