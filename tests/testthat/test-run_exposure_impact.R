test_that("run_exposure_impact computes exposure-level centrality scores", {
    library(MultiAssayExperiment)

    # Create dummy exposomicset
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # mae differential abundance
    mae <- run_differential_abundance(
        exposomicset = mae,
        formula = ~ smoker + sex,
        abundance_col = "counts",
        method = "limma_voom",
        action = "add"
    )

    # Run DEG association + correlation
    mae <- run_correlation(mae,
        feature_type = "degs",
        exposure_cols = c("age", "bmi", "exposure_pm25"),
        correlation_cutoff = 0,
        feature_cors = T,
        deg_pval_thresh = 0.05,
        deg_pval_col = "P.Value",
        pval_cutoff = 0.1,
        action = "add"
    )

    mae <- run_correlation(mae,
        feature_type = "degs",
        exposure_cols = c("age", "bmi", "exposure_pm25"),
        correlation_cutoff = 0,
        feature_cors = F,
        deg_pval_thresh = 0.05,
        deg_pval_col = "P.Value",
        pval_cutoff = 0.1,
        action = "add"
    )

    # Build networks
    mae <- run_create_network(mae,
        feature_type = "degs",
        action = "add"
    )

    mae <- run_create_network(mae,
        feature_type = "degs_feature_cor",
        action = "add"
    )

    # Run exposure impact
    mae <- run_exposure_impact(mae,
        feature_type = "degs",
        action = "add"
    )

    # Extract result
    impact <- metadata(mae)$network$exposure_impact$degs$exposure_impact

    # Expectations
    # ensure it is a data frame
    expect_s3_class(impact, "data.frame")

    # ensure the expected columns are present
    expect_true(all(c("exposure", "mean_degree", "n_features") %in% colnames(impact)))
})
