test_that("filter_non_normal removes non-normal exposures in expomicset", {
  # Create dummy data
  dummy <- make_dummy_data(n_samples = 100)

  # Create MultiAssayExperiment
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # run the normality check
  mae <- mae |>
    run_normality_check() |>
    transform_exposure()

  # Ensure normality metadata is present
  norm_df <- metadata(mae)$quality_control$transformation$norm_df
  expect_true(is.data.frame(norm_df))
  expect_true(all(c("exposure", "p.value") %in% colnames(norm_df)))

  # Identify non-normal exposures
  non_normal <- norm_df$exposure[norm_df$p.value < 0.05]

  # Apply filter
  filtered <- filter_non_normal(mae, p_thresh = 0.05)

  # Ensure non-normal exposures are removed from main colData
  expect_false(any(non_normal %in% colnames(colData(filtered))))

  # Ensure non-normal exposures are removed from all experiment colData
  for (exp_name in names(experiments(filtered))) {
    expect_false(any(non_normal %in% colnames(colData(experiments(filtered)[[exp_name]]))))
  }

  # Confirm metadata was updated
  steps <- metadata(filtered)$summary$steps
  expect_true("filter_non_normal" %in% names(steps))
  expect_match(steps$filter_non_normal$notes, "Filtered out non-normal exposure variables")
})
