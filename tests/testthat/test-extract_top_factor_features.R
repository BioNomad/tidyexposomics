test_that("extract_top_factor_features works for MOFA", {
    # skip unless MOFA2 is available and Python backend is good to go
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        skip("MOFA2 not available")
    }
    # skip on arm64 darwin macs where MOFA2 breaks
    if (Sys.info()[["sysname"]] == "Darwin" &&
        grepl("arm64", R.version$platform)) {
        skip("Skip MOFA2 test on macOS ARM64 - issue with CI segfault")
    }

    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # with this small of a dataset it needs to be 2 factors
    mae <- run_multiomics_integration(mae, method = "MOFA", n_factors = 2)

    factors <- colnames(MOFA2::get_factors(metadata(mae)$multiomics_integration$integration_results$result)[[1]])

    top_feats <- extract_top_factor_features(
        mae,
        factors = factors,
        method = "percentile",
        percentile = 0.9,
        action = "get"
    )

    # does this return a data frame
    expect_s3_class(top_feats, "data.frame")

    # do we get the expected column names
    expect_true(all(c("feature", "factor", "loading", "exp_name") %in% colnames(top_feats)))

    # do we have more than 0 rows for the top features
    expect_gt(nrow(top_feats), 0)
})

test_that("extract_top_factor_features works for MCIA", {
    skip_if_not_installed("nipalsMCIA")

    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    mae <- run_multiomics_integration(mae, method = "MCIA", n_factors = 3)

    factors <- mae@metadata$multiomics_integration$integration_results$result@global_scores |>
        as.data.frame() |>
        colnames()

    top_feats <- extract_top_factor_features(
        mae,
        factors = factors,
        method = "percentile",
        percentile = 0.9,
        action = "get"
    )

    # does this return a data frame
    expect_s3_class(top_feats, "data.frame")

    # do we get the expected column names
    expect_true(all(c("feature", "factor", "loading", "exp_name") %in% colnames(top_feats)))

    # do we have more than 0 rows for the top features
    expect_gt(nrow(top_feats), 0)
})

test_that("extract_top_factor_features works for DIABLO", {
    skip_if_not_installed("mixOmics")

    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    mae <- run_multiomics_integration(
        mae,
        method = "DIABLO",
        outcome = "smoker",
        n_factors = 3
    )

    # adding the assay name for multi-block components
    factors <- colnames(metadata(mae)$multiomics_integration$integration_results$result$variates[["mRNA"]])
    factors <- paste0("mRNA ", factors)

    top_feats <- extract_top_factor_features(
        mae,
        factors = factors,
        method = "percentile",
        percentile = 0.9,
        action = "get"
    )

    # does this return a data frame
    expect_s3_class(top_feats, "data.frame")

    # do we get the expected column names
    expect_true(all(c("feature", "factor", "loading", "exp_name") %in% colnames(top_feats)))

    # do we have more than 0 rows for the top features
    expect_gt(nrow(top_feats), 0)
})

test_that("extract_top_factor_features works for RGCCA", {
    skip_if_not_installed("RGCCA")

    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    mae <- run_multiomics_integration(mae,
        method = "RGCCA",
        n_factors = 3
    )

    # adding the assay name for multi-block components
    factors <- colnames(metadata(mae)$multiomics_integration$integration_results$result$a[["mRNA"]] |>
        as.data.frame())
    factors <- paste0("mRNA ", factors)

    top_feats <- extract_top_factor_features(
        mae,
        factors = factors,
        method = "percentile",
        percentile = 0.9,
        action = "get"
    )

    # does this return a data frame
    expect_s3_class(top_feats, "data.frame")

    # do we get the expected column names
    expect_true(all(c("feature", "factor", "loading", "exp_name") %in% colnames(top_feats)))

    # do we have more than 0 rows for the top features
    expect_gt(nrow(top_feats), 0)
})
