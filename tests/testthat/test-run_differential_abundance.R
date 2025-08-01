test_that("run_differential_abundance returns expected results and updates metadata", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 20)

    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Run differential abundance
    mae_da <- run_differential_abundance(
        expomicset = mae,
        formula = ~ smoker + sex,
        abundance_col = "counts",
        method = "limma_voom",
        action = "add"
    )

    # Check that results were added to metadata
    da_results <- MultiAssayExperiment::metadata(mae_da)$differential_analysis$differential_abundance

    # is it the right class
    expect_s3_class(da_results, "tbl_df")

    # do we get the expected column names
    expect_true(all(c("feature", "logFC", "adj.P.Val", "exp_name") %in% colnames(da_results)))

    # Confirm a step was recorded
    steps <- MultiAssayExperiment::metadata(mae_da)$summary$steps
    expect_true("run_differential_abundance" %in% names(steps))
    expect_match(steps$run_differential_abundance$notes, "Performed differential abundance analysis")
})
