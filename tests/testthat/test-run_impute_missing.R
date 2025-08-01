test_that("run_impute_missing() performs exposure and omics imputation correctly", {
    skip_if_not_installed("naniar")
    skip_if_not_installed("impute")

    # Create dummy data with missing values injected
    dummy <- make_example_data(n_samples = 20)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Inject NA values into exposure and omics
    colData(mae)$exposure_pm25[1:3] <- NA
    assay(mae[["mRNA"]])[1:2, 1:2] <- NA

    # Add minimal NA QC metadata to trigger exposure imputation
    metadata(mae)$quality_control <- list(
        na_qc = list(exposure = TRUE)
    )

    # Run imputation
    imputed <- run_impute_missing(
        mae,
        exposure_impute_method = "median",
        exposure_cols = c("exposure_pm25"),
        omics_impute_method = "knn",
        omics_to_impute = c("mRNA")
    )

    # Exposure imputation, no NAs should remain in selected column
    expect_false(any(is.na(colData(imputed)$exposure_pm25)))

    # Omics imputation, no NAs should remain in imputed assay
    imputed_assay <- assay(experiments(imputed)[["mRNA"]])
    expect_false(any(is.na(imputed_assay)))

    # Step metadata should be updated
    step_log <- metadata(imputed)$summary$steps
    expect_true(any(names(step_log) == "run_impute_missing"))
})
