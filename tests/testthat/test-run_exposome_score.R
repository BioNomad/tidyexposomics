test_that("run_exposome_score computes scores with all supported methods", {
  skip_if_not_installed("matrixStats")

  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # try all the methods
  methods <- c("mean", "sum", "median", "pca", "quantile", "var")

  for (method in methods) {
    result <- run_exposome_score(
      expomicset = mae,
      score_type = method,
      scale = TRUE
    )
    # what is the score column
    score_col <- grep("exposome_score", colnames(colData(result)), value = TRUE)

    # is the score column there?
    expect_true(length(score_col) == 1, label = paste("Missing score column for", method))

    # esnure that it is the right type
    expect_type(colData(result)[[score_col]], "double")

    # make sure it is equal to the number of samples
    expect_equal(nrow(colData(result)), 20)
  }
})

test_that("run_exposome_score computes score for user-defined columns", {
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  exposure_vars <- c("age","bmi","exposure_pm25")

  result <- run_exposome_score(
    expomicset = mae,
    score_type = "sum",
    exposure_cols = exposure_vars,
    scale = FALSE,
    score_column_name = "custom_score"
  )

  # ensure the column name is there
  expect_true("custom_score" %in% colnames(colData(result)))

  # ensure it is the right type
  expect_type(colData(result)$custom_score, "double")
})

test_that("run_exposome_score with IRT works if mirt is available", {
  skip_if_not_installed("mirt")

  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  exposure_vars <- c("age","bmi","exposure_pm25")

  result <- run_exposome_score(
    expomicset = mae,
    score_type = "irt",
    exposure_cols = exposure_vars,
    scale = TRUE
  )

  # ensure column is present
  expect_true("exposome_score_irt" %in% colnames(colData(result)))

  # ensure column is the right type
  expect_type(colData(result)$exposome_score_irt, "double")
})

test_that("run_exposome_score throws error for invalid method", {
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # ensure it throws an error with an unsupported score_type
  expect_error(
    run_exposome_score(mae, score_type = "invalid"),
    "Invalid score_type"
  )
})
