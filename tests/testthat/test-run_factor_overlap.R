test_that("run_factor_overlap() returns expected object and annotations", {
    skip_if_not_installed("MultiAssayExperiment")

    # Create dummy data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Run differential abundance
    mae <- run_differential_abundance(
        exposomicset = mae,
        formula = ~ smoker + sex,
        abundance_col = "counts",
        method = "limma_voom",
        action = "add"
    )

    # Perform multi-omics integration and extract top features
    mae <- run_multiomics_integration(mae,
        method = "DIABLO",
        outcome = "smoker",
        n_factors = 3
    )

    mae <- extract_top_factor_features(mae,
        percentile = .5,
        factors = c("V1", "V2", "V3")
    )

    # Run function
    mae <- run_factor_overlap(mae,
        pval_col = "P.Value",
        robust = F,
        action = "add"
    )

    # Check output structure
    expect_s4_class(mae, "MultiAssayExperiment")

    # are the common factor features added to the metadata
    expect_true("common_top_factor_features" %in%
        names(mae@metadata$multiomics_integration))

    res_df <- run_factor_overlap(mae,
        robust = F,
        pval_col = "P.Value",
        action = "get"
    )

    # is the result a data frame
    expect_s3_class(res_df, "data.frame")

    # are the expected column names present
    expect_true(all(c("exp_name", "feature", "factor", "is_deg") %in% names(res_df)))

    # check if the is_deg column populated correctly
    expect_true(is.logical(res_df$is_deg))
})
