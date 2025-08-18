test_that("filter_missing removes features and generates QC plots", {
    skip_if_not_installed("naniar")

    # Create dummy data with injected missingness
    dummy <- make_example_data(n_samples = 10)
    dummy$exposure$age[1:3] <- NA
    dummy$omics$mRNA[1:2, 1:5] <- NA
    dummy$omics$proteomics[1:4, 1:3] <- NA

    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # filter it so that variables with more than 20% NAs are filtered out
    mae_filtered <- filter_missing(mae, na_thresh = 20)

    # Test filtered exposure columns
    # get the missing summary for the coldata
    col_missing <- naniar::miss_var_summary(as.data.frame(colData(mae)))

    # what do we expect to be removed
    expected_removed <- col_missing$variable[col_missing$pct_miss > 20]

    # we want this to be false becuase the removed variables should not be
    # in the coldata of the filtered object
    expect_false(any(expected_removed %in% colnames(colData(mae_filtered))))

    # Test filtered omics features
    for (omic in names(experiments(mae))) {
        # check the rownames of the original data and filtered
        raw_feats <- rownames(dummy$omics[[omic]])
        filtered_feats <- rownames(mae_filtered[[omic]])

        # transpose the omics data
        assay_data <- t(dummy$omics[[omic]])

        # generate the missinness summary
        miss_summary <- naniar::miss_var_summary(as.data.frame(assay_data))

        # what do we expect to be removed
        expected_removed <- miss_summary$variable[miss_summary$pct_miss > 20]

        # ensure that the features with too many NAs are out of the object
        expect_false(any(expected_removed %in% filtered_feats))
    }

    # Check that QC plots exist in metadata
    qc_meta <- metadata(mae_filtered)$quality_control$na_qc
    expect_type(qc_meta$exposure$na_exposure_qc_plot, "list")
    for (omic in names(experiments(mae_filtered))) {
        expect_type(qc_meta[[omic]]$na_omics_qc_plot, "list")
    }
})
