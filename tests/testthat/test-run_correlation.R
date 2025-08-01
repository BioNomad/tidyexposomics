test_that("run_correlation computes correlations between exposures and omics", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 30)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Run correlation
    correlated <- run_correlation(
        expomicset = mae,
        feature_type = "omics",
        correlation_method = "spearman",
        correlation_cutoff = 0.1,
        pval_cutoff = 1,
        action = "add"
    )

    # Confirm metadata contains the correlation results
    corr_meta <- MultiAssayExperiment::metadata(correlated)$correlation$omics

    # does it deliver a data frame
    expect_s3_class(corr_meta, "data.frame")

    # are all the expected columns present
    expect_true(all(c("exposure", "feature", "correlation", "p.value") %in% colnames(corr_meta)))

    # Confirm a record was added to the summary steps
    steps <- MultiAssayExperiment::metadata(correlated)$summary$steps

    # is the step added
    expect_true("run_correlation_omics" %in% names(steps))
    expect_match(steps$run_correlation_omics$notes, "Correlated omics features with exposures")
})

test_that("run_correlation computes correlations between exposures", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 30)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Run correlation
    correlated <- run_correlation(
        expomicset = mae,
        feature_type = "exposures",
        correlation_method = "spearman",
        correlation_cutoff = 0.1,
        pval_cutoff = 1,
        action = "add"
    )

    # Confirm metadata contains the correlation results
    corr_meta <- MultiAssayExperiment::metadata(correlated)$correlation$exposures

    # does it deliver a data frame
    expect_s3_class(corr_meta, "data.frame")

    # are all the expected columns present
    expect_true(all(c("var1", "var2", "correlation", "p.value") %in% colnames(corr_meta)))

    # Confirm a record was added to the summary steps
    steps <- MultiAssayExperiment::metadata(correlated)$summary$steps

    # is the step added
    expect_true("run_correlation_exposures" %in% names(steps))
    expect_match(steps$run_correlation_exposures$notes, "Correlated exposures features with exposures")
})
