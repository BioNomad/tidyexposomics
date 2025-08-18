test_that("run_cluster_samples performs clustering and records output", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # Choose two numeric exposures to cluster on
    available_numeric <- mae |>
        MultiAssayExperiment::colData() |>
        as.data.frame() |>
        dplyr::select(where(is.numeric)) |>
        colnames()

    expect_true(length(available_numeric) >= 2)

    selected_exposures <- available_numeric[1:2]

    # Run clustering
    clustered <- run_cluster_samples(
        expomicset = mae,
        exposure_cols = selected_exposures,
        dist_method = "euclidean",
        clustering_approach = "diana",
        cluster_method = "ward.D",
        action = "add"
    )

    # Check that sample_group column is added to colData
    expect_true("sample_group" %in% colnames(MultiAssayExperiment::colData(clustered)))

    # Check that metadata contains clustering output
    qc_meta <- MultiAssayExperiment::metadata(clustered)$quality_control$sample_clustering
    expect_true("sample_cluster" %in% names(qc_meta))
    expect_true("sample_groups" %in% names(qc_meta))
    expect_s3_class(qc_meta$sample_cluster, "hclust")

    # Check that step was recorded
    steps <- MultiAssayExperiment::metadata(clustered)$summary$steps
    expect_true("run_cluster_samples" %in% names(steps))
    expect_match(steps$run_cluster_samples$notes, "Optimal number of clusters")
})
