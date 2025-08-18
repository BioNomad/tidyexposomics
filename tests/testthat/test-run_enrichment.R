test_that("run_enrichment identifies GO terms for DEGs and updates metadata", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 30)

    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Run differential abundance to generate DEGs
    mae <- run_differential_abundance(
        expomicset = mae,
        formula = ~ smoker + sex,
        abundance_col = "counts",
        method = "limma_voom",
        action = "add"
    )

    # Run enrichment on DEGs using GO
    enriched <- run_enrichment(
        expomicset = mae,
        feature_type = "degs",
        feature_col = "symbol",
        deg_pval_col = "P.Value",
        deg_pval_threshold = 0.5,
        clustering_approach = "diana",
        deg_logfc_col = "logFC",
        deg_logfc_threshold = log2(1.5),
        db = "GO",
        species = "goa_human",
        action = "add"
    )

    # Confirm enrichment results added to metadata
    enr_res <- MultiAssayExperiment::metadata(enriched)$enrichment$degs
    expect_s3_class(enr_res, "data.frame")
    expect_true(all(c("term_id", "term_name", "p_value", "padj", "n_with_sel", "ids") %in% colnames(enr_res)))


    # Confirm a step record was logged
    steps <- MultiAssayExperiment::metadata(enriched)$summary$steps
    expect_true("run_enrichment" %in% names(steps))
    expect_match(steps$run_enrichment$notes, "Performed GO enrichment on degs features")
})
