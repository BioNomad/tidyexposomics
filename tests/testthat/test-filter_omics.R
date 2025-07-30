test_that("filter_omics removes low-quality features using variance method", {
  # Create dummy data
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Get initial feature counts
  initial_nrow <- nrow(assay(experiments(mae)[["proteomics"]]))

  # Artificially lower expression in proteomics
  assay(experiments(mae)[["proteomics"]])[1:15, ] <- 0

  # Apply variance filter
  filtered <- filter_omics(
    expomicset = mae,
    method = "variance",
    assays = "proteomics",
    assay_name = 1,
    min_var = 0.8,
    verbose = FALSE
  )

  # Check that features were removed
  final_nrow <- nrow(assay(experiments(filtered)[["proteomics"]]))
  expect_lt(final_nrow, initial_nrow)

  # Check metadata logging
  step_names <- names(metadata(filtered)$summary$steps)
  expect_true(any(grepl("filter_omics_proteomics", step_names)))
})

test_that("filter_omics removes low-expressed features using expression method", {
  # Create dummy data
  dummy <- make_example_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Artificially lower expression in mRNA
  assay(experiments(mae)[["mRNA"]])[1:5, ] <- 0

  # Get initial feature counts
  initial_nrow <- nrow(assay(experiments(mae)[["mRNA"]]))

  # Apply expression filter
  filtered <- filter_omics(
    expomicset = mae,
    method = "expression",
    assays = "mRNA",
    assay_name = 1,
    min_value = 1,
    min_prop = 0.5,
    verbose = FALSE
  )

  # Check that features were removed
  final_nrow <- nrow(assay(experiments(filtered)[["mRNA"]]))
  expect_lt(final_nrow, initial_nrow)

  # Check metadata logging
  step_names <- names(metadata(filtered)$summary$steps)
  expect_true(any(grepl("filter_omics_mRNA", step_names)))
})
