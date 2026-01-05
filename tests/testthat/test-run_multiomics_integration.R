test_that("run_multiomics_integration works with MOFA", {
    # skipping on Bioconductor because of basilisk issues
    skip_on_bioc()
    # MOFA2 has an issue running on Mac-arm64
    skip_on_os("mac")
    skip_if(grepl("arm64", R.version$platform) &&
        Sys.info()[["sysname"]] == "Darwin")
    skip_if_not_installed("MOFA2")

    local({
        dummy <- make_example_data(n_samples = 20)
        mae <- create_exposomicset(
            codebook = dummy$codebook,
            exposure = dummy$exposure,
            omics = dummy$omics,
            row_data = dummy$row_data
        )

        mae <- run_multiomics_integration(
            mae,
            method = "MOFA",
            n_factors = 3,
            scale = TRUE,
            action = "add"
        )

        res <- metadata(mae)$multiomics_integration$integration_results

        # check that the right method was recorded
        expect_equal(res$method, "MOFA")

        # safely grab MOFA2 functions
        get_factors <- get("get_factors", asNamespace("MOFA2"))
        get_weights <- get("get_weights", asNamespace("MOFA2"))

        # ensure that sample level scores are captured
        expect_equal(
            rownames(get_factors(res$result)[[1]]) |> sort(),
            rownames(mae@colData) |> sort()
        )

        # ensure that loading scores are captured
        for (exp_name in names(mae@ExperimentList)) {
            expect_equal(
                rownames(mae@ExperimentList[[exp_name]]),
                rownames(get_weights(res$result)[[exp_name]])
            )
        }

        # ensure that the step was recorded
        step_names <- names(metadata(mae)$summary$steps)
        expect_true(any(grepl("run_multiomics_integration", step_names)))
    })
})


test_that("run_multiomics_integration works with MCIA", {
    # Skipping on Bioc due to RMD check failure
    skip_on_bioc()
    skip_if_not_installed("nipalsMCIA")

    # make the example data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run multiomics integration
    mae <- run_multiomics_integration(
        mae,
        method = "MCIA",
        n_factors = 3,
        scale = TRUE,
        action = "add"
    )

    # grab the integration results
    res <- metadata(mae)$multiomics_integration$integration_results

    # check that the right method was recorded
    expect_equal(res$method, "MCIA")

    # ensure that sample level scores are captured
    # by seeing if the sample ids are present
    expect_equal(
        rownames(res$result@global_scores) |>
            sort(),
        rownames(mae@colData) |>
            sort()
    )

    # ensure that loading scores are captured
    for (exp_name in names(mae@ExperimentList)) {
        expect_equal(
            rownames(mae@ExperimentList[[exp_name]]),
            rownames(res$result@block_loadings[[exp_name]])
        )
    }

    # ensure that the step was recorded
    step_names <- names(metadata(mae)$summary$steps)
    expect_true(any(grepl("run_multiomics_integration", step_names)))
})


test_that("run_multiomics_integration works with RGCCA", {
    skip_if_not_installed("RGCCA")

    # make the example data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run multiomics integration
    mae <- run_multiomics_integration(
        mae,
        method = "RGCCA",
        n_factors = 3,
        scale = TRUE,
        action = "add"
    )

    # grab the integration results
    res <- metadata(mae)$multiomics_integration$integration_results

    # check that the right method was recorded
    expect_equal(res$method, "RGCCA")


    # ensure that sample and loading scores are captured
    for (exp_name in names(mae@ExperimentList)) {
        # ensure that sample level scores are captured
        # by seeing if the sample ids are present
        # doing this inside the loop since these are calculated
        # on a block-wise basis
        expect_equal(
            rownames(res$result$Y[[exp_name]]) |>
                sort(),
            rownames(mae@colData) |>
                sort()
        )
        # ensure loading scores are present
        expect_equal(
            rownames(mae@ExperimentList[[exp_name]]),
            rownames(res$result$a[[exp_name]])
        )
    }

    # ensure that the step was recorded
    step_names <- names(metadata(mae)$summary$steps)
    expect_true(any(grepl("run_multiomics_integration", step_names)))
})


test_that("run_multiomics_integration works with DIABLO", {
    skip_if_not_installed("mixOmics")

    # make the example data
    dummy <- make_example_data(n_samples = 20)
    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # run multiomics integration
    mae <- run_multiomics_integration(
        mae,
        method = "DIABLO",
        outcome = "smoker",
        n_factors = 3,
        scale = TRUE,
        action = "add"
    )

    # grab the integration results
    res <- metadata(mae)$multiomics_integration$integration_results

    # check that the right method was recorded
    expect_equal(res$method, "DIABLO")


    # ensure that sample and loading scores are captured
    for (exp_name in names(mae@ExperimentList)) {
        # ensure that sample level scores are captured
        # by seeing if the sample ids are present
        # doing this inside the loop since these are calculated
        # on a block-wise basis
        expect_equal(
            rownames(res$result$variates[[exp_name]]) |>
                sort(),
            rownames(mae@colData) |>
                sort()
        )
        # ensure loading scores are present
        expect_equal(
            rownames(mae@ExperimentList[[exp_name]]),
            rownames(res$result$loadings[[exp_name]])
        )
    }

    # ensure that the step was recorded
    step_names <- names(metadata(mae)$summary$steps)
    expect_true(any(grepl("run_multiomics_integration", step_names)))
})
