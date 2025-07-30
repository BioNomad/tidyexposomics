test_that("run_association works with exposures as source and guassian", {
  # Create dummy data
  dummy <- make_example_data(n_samples = 30)

  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Select a few exposures to associate
  feature_set <- intersect(c("age", "bmi", "alcohol"), colnames(dummy$exposure))

  # Run association analysis
  assoc_mae <- run_association(
    expomicset = mae,
    outcome = "exposure_pm25",
    source = "exposures",
    covariates = "sex",
    feature_set = feature_set,
    family = "gaussian",
    action = "add"
  )

  # Confirm metadata was added
  assoc_meta <- MultiAssayExperiment::metadata(assoc_mae)$association$assoc_exposures
  expect_true("results_df" %in% names(assoc_meta))
  expect_true("covariates" %in% names(assoc_meta))
  expect_true(all(c("term", "estimate", "p.value", "r2", "adj_r2") %in% colnames(assoc_meta$results_df)))

  # Confirm step was recorded
  steps <- MultiAssayExperiment::metadata(assoc_mae)$summary$steps
  expect_true("run_association" %in% names(steps))
  expect_match(steps$run_association$notes, "Performed association analysis")
})


test_that("run_association works with exposures as source and binomial", {
  # Create dummy data
  dummy <- make_example_data(n_samples = 30)

  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Select a few exposures to associate
  feature_set <- intersect(c("age", "bmi", "alcohol"), colnames(dummy$exposure))

  # Run association analysis
  assoc_mae <- run_association(
    expomicset = mae,
    outcome = "smoker",
    source = "exposures",
    covariates = "sex",
    feature_set = feature_set,
    family = "binomial",
    action = "add"
  )

  # Confirm metadata was added
  assoc_meta <- MultiAssayExperiment::metadata(assoc_mae)$association$assoc_exposures
  expect_true("results_df" %in% names(assoc_meta))
  expect_true("covariates" %in% names(assoc_meta))
  expect_true(all(c("term", "estimate", "p.value", "r2", "adj_r2") %in% colnames(assoc_meta$results_df)))

  # Confirm step was recorded
  steps <- MultiAssayExperiment::metadata(assoc_mae)$summary$steps
  expect_true("run_association" %in% names(steps))
  expect_match(steps$run_association$notes, "Performed association analysis")
})

