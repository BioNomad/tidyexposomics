test_that("extract_omics_exposure_df works with and without variable_map", {
    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Without variable_map
    merged_df <- extract_omics_exposure_df(mae)

    # is it a df
    expect_s3_class(merged_df, "data.frame")

    # are the rownames of the returned df the same as the rownames of the
    # original meta data
    expect_true(all(rownames(merged_df) %in% rownames(dummy$exposure)))

    # ensure that the colnames have the right prefixes
    # should be features
    expect_true(any(grepl("mRNA_", colnames(merged_df))))
    expect_true(any(grepl("proteomics_", colnames(merged_df))))

    # ensure the number of columns is greater
    # than the number of exposures alone
    expect_gt(ncol(merged_df), 5)

    # With variable_map
    gene_subset <- rownames(SummarizedExperiment::assay(mae[["mRNA"]]))[1:5]
    prot_subset <- rownames(SummarizedExperiment::assay(mae[["proteomics"]]))[1:5]
    var_map <- rbind(
        data.frame(variable = gene_subset, exp_name = "mRNA"),
        data.frame(variable = prot_subset, exp_name = "proteomics"),
        data.frame(variable = c("age", "bmi"), exp_name = "exposure")
    )

    merged_df_map <- extract_omics_exposure_df(mae, variable_map = var_map)

    # ensure we have the right prefixes - again features
    expect_true(any(grepl("^mRNA_", colnames(merged_df_map))))
    expect_true(any(grepl("^proteomics_", colnames(merged_df_map))))

    # make sure that we have our metadata variables in the returned
    # df
    expect_true(all(c("age", "bmi") %in% colnames(merged_df_map)))
})

test_that("extract_omics_exposure_df handles log2 transformation correctly", {
    dummy <- make_example_data(n_samples = 10)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    df_log2 <- extract_omics_exposure_df(mae, log2_trans = TRUE)
    df_raw <- extract_omics_exposure_df(mae, log2_trans = FALSE)

    m_col <- grep("^mRNA_", colnames(df_log2), value = TRUE)[1]

    # ensure the log2 transformed values are less than the raw values
    expect_lt(mean(df_log2[[m_col]]), mean(df_raw[[m_col]]))
})
