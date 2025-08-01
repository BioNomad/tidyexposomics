test_that("run_normality_check works as expected", {
    # make the example data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run the normality check
    mae_checked <- run_normality_check(mae,
        action = "add"
    )

    # grab the results from metadata
    norm_meta <- metadata(mae_checked)$quality_control$normality

    # check that the norm_df exists
    expect_s3_class(norm_meta$norm_df, "data.frame")

    # are the expected column names present
    expect_true(all(c("statistic", "p.value", "method", "exposure") %in% colnames(norm_meta$norm_df)))

    # check that the summary is a data frame with correct structure
    expect_s3_class(norm_meta$norm_summary, "data.frame")
    expect_true(all(c("var", "value") %in% colnames(norm_meta$norm_summary)))

    # check that the step was recorded
    step_names <- names(metadata(mae_checked)$summary$steps)
    expect_true("run_normality_check" %in% step_names)
})
