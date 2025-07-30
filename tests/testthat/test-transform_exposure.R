test_that("transform_exposure applies boxcox_best correctly", {
  skip_if_not_installed("MASS")

  # Create dummy data
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Select a few exposures to transform
  exposures <- colnames(SummarizedExperiment::colData(mae))[1:3]

  # Apply transformation
  transformed <- transform_exposure(
    expomicset = mae,
    exposure_cols = exposures,
    transform_method = "boxcox_best"
  )

  # Check that colData still contains the transformed columns
  transformed_data <- as.data.frame(SummarizedExperiment::colData(transformed))
  expect_true(all(exposures %in% colnames(transformed_data)))

  # Check that metadata has transformation results
  trans_meta <- MultiAssayExperiment::metadata(transformed)$quality_control$transformation
  expect_s3_class(trans_meta$norm_df, "data.frame")
  expect_true("p.value" %in% colnames(trans_meta$norm_df))
  expect_true("exposure" %in% colnames(trans_meta$norm_df))

  # Check codebook is updated
  codebook <- MultiAssayExperiment::metadata(transformed)$codebook
  expect_s3_class(codebook, "data.frame")
  expect_true(all(exposures %in% codebook$variable))
  expect_true("transformation_applied" %in% colnames(codebook))

  # Check step summary
  steps <- MultiAssayExperiment::metadata(transformed)$summary$steps
  expect_true("transform_exposure" %in% names(steps))
  expect_match(steps$transform_exposure$notes, "Applied 'boxcox_best' transformation")
})


test_that("transform_exposure applies log2 transformation", {
  # Create dummy data
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Select a few exposures to transform
  exposures <- c("age","bmi","exposure_pm25")

  # Apply transformation
  transformed <- transform_exposure(
    expomicset = mae,
    exposure_cols = exposures,
    transform_method = "log2"
  )

  # Check that colData still contains the transformed columns
  transformed_data <- as.data.frame(SummarizedExperiment::colData(transformed))
  expect_true(all(exposures %in% colnames(transformed_data)))

  # Check transformation metadata is present
  trans_meta <- MultiAssayExperiment::metadata(transformed)$quality_control$transformation
  expect_s3_class(trans_meta$norm_df, "data.frame")
  expect_true("p.value" %in% colnames(trans_meta$norm_df))
  expect_true("exposure" %in% colnames(trans_meta$norm_df))

  # Check transformation label is "log2"
  codebook <- MultiAssayExperiment::metadata(transformed)$codebook
  expect_true(all(codebook[codebook$variable %in% exposures,]$transformation_applied == "log2"))

  # Check step summary
  steps <- MultiAssayExperiment::metadata(transformed)$summary$steps
  expect_true("transform_exposure" %in% names(steps))
  expect_match(steps$transform_exposure$notes, "Applied 'log2' transformation")
})
