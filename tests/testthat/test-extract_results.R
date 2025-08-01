test_that("extract_results correctly retrieves results from metadata", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 5)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Inject mock results into metadata
    metadata(mae)$association$omics <- data.frame(
        feature = "B2M",
        exposure = "age",
        p = 0.01
    )

    metadata(mae)$enrichment$omics <- data.frame(
        term = "antigen presentation",
        p_value = 0.01,
        ids = "HLA-DRB1,TAP1,B2M"
    )

    # Extract and check results
    assoc_result <- extract_results(mae, result = "association")
    enrich_result <- extract_results(mae, result = "enrichment")

    # ensure it is pulling a df and that the expected column names are present
    expect_s3_class(assoc_result$omics, "data.frame")
    expect_named(assoc_result$omics, c("feature", "exposure", "p"))
    expect_s3_class(enrich_result$omics, "data.frame")
    expect_true("term" %in% names(enrich_result$omics))
})



test_that("extract_results returns NULL for missing metadata entry", {
    dummy <- make_example_data(n_samples = 5)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Should return NULL for result not present
    missing_result <- extract_results(mae, result = "network")
    expect_null(missing_result)
})
