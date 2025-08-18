test_that("filter_sample_outliers removes PCA-identified outliers", {
    # make the example data
    dummy <- make_example_data(n_samples = 20)

    # artificially make a sample an outlier
    dummy$omics$mRNA[1:3, 1] <- c(1e5, -1e5, 1e5)

    # create the multiassayexperiment
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run PCA to detect outliers
    mae_pca <- run_pca(mae,
        action = "add"
    )

    # get initial sample count
    initial_n <- nrow(colData(mae_pca))

    # run outlier filter
    mae_filtered <- filter_sample_outliers(mae_pca)

    # check that sample count decreased (or stayed the same)
    expect_lt(nrow(colData(mae_filtered)), initial_n)

    # check that PCA outliers are no longer in colData
    outliers <- metadata(mae_pca)$quality_control$pca$outliers
    expect_false(any(outliers %in% rownames(colData(mae_filtered))))

    # ensure that the step was recorded
    step_names <- names(metadata(mae_filtered)$summary$steps)
    expect_true("filter_sample_outliers" %in% step_names)
})
