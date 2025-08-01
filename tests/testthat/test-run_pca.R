test_that("run_pca works as expected", {
    # make the example data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run PCA
    mae_pca <- run_pca(mae,
        log_trans_exp = TRUE,
        log_trans_omics = TRUE,
        action = "add"
    )

    # grab the PCA results from metadata
    pca_meta <- metadata(mae_pca)$quality_control$pca

    # check that PCA results exist and are structured correctly
    expect_s3_class(pca_meta$pca_df, "tbl_df")
    expect_s3_class(pca_meta$pca_sample, "prcomp")
    expect_s3_class(pca_meta$pca_feature, "prcomp")
    expect_type(pca_meta$outliers, "character")

    # grab the colData
    col_data <- colData(mae_pca) |>
        as.data.frame()

    # ensure that PCs are added to colData
    expect_true(any(grepl("^PC[0-9]+$", colnames(col_data))))

    # ensure that the step was recorded
    step_names <- names(metadata(mae_pca)$summary$steps)
    expect_true("run_pca" %in% step_names)
})
