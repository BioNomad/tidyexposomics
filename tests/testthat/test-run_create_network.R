test_that("run_create_network builds igraph object from omics_feature_cor", {
    # Create dummy data
    dummy <- make_example_data(n_samples = 30)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # compute feature-feature correlations (required for network)
    mae <- run_correlation(
        expomicset = mae,
        feature_type = "omics",
        feature_cors = TRUE,
        correlation_method = "spearman",
        correlation_cutoff = 0.1,
        pval_cutoff = 1,
        action = "add"
    )

    # build correlation network
    mae_net <- run_create_network(
        expomicset = mae,
        feature_type = "omics_feature_cor",
        action = "add"
    )

    # Check that network was added to metadata
    net_obj <- MultiAssayExperiment::metadata(mae_net)$network$network_omics_feature_cor

    # make sure the graph and graph summary stats are there
    expect_true("graph" %in% names(net_obj))
    expect_true("summary" %in% names(net_obj))

    # make sure the graph and summary stats are the right class
    expect_s3_class(net_obj$graph, "igraph")
    expect_s3_class(net_obj$summary, "tbl_df")

    # Check that analysis step was recorded
    steps <- MultiAssayExperiment::metadata(mae_net)$summary$steps
    expect_true("run_create_network" %in% names(steps))
    expect_match(steps$run_create_network$notes, "Created undirected network")
})

test_that("run_create_network returns igraph object when action = 'get'", {
    dummy <- make_example_data(n_samples = 30)
    mae <- create_expomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    mae <- run_correlation(
        expomicset = mae,
        feature_type = "omics",
        feature_cors = TRUE,
        correlation_method = "spearman",
        correlation_cutoff = 0.1,
        pval_cutoff = 1,
        action = "add"
    )

    net <- run_create_network(
        expomicset = mae,
        feature_type = "omics_feature_cor",
        action = "get"
    )

    # ensure that we get the graph returned
    expect_true("graph" %in% names(net))

    # ensure the class is correct
    expect_s3_class(net$graph, "igraph")
})
