# Load Libraries
library(testthat)
library(MultiAssayExperiment)
library(S4Vectors)

test_that("create_exposomicset returns a valid MultiAssayExperiment", {
    dummy <- make_example_data(n_samples = 10, n_proteins = 50)

    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics,
        row_data = dummy$row_data
    )

    # is it a MAE
    expect_s4_class(mae, "MultiAssayExperiment")

    # do we have the expected assays
    expect_named(experiments(mae), c("mRNA", "proteomics"))

    # do we have the same samples in coldata and exposure
    expect_equal(nrow(colData(mae)), nrow(dummy$exposure))

    # do the features line up
    expect_equal(nrow(mae[["mRNA"]]), nrow(dummy$omics$mRNA))
    expect_equal(nrow(mae[["proteomics"]]), nrow(dummy$omics$proteomics))

    # do the rownames line up
    expect_equal(rownames(rowData(mae[["mRNA"]])), rownames(dummy$omics$mRNA))

    # does the feature data line up
    expect_equal(rowData(mae[["mRNA"]])$symbol[1], dummy$row_data$mRNA$symbol[1])

    # do we have the same samples as in the exposure df
    expect_true(all((colnames(mae) |>
        unlist() |>
        unique()) %in% rownames(dummy$exposure)))

    # did we add the codebook
    expect_true("codebook" %in% names(metadata(mae)))
})

test_that("create_exposomicset works with single omics input", {
    dummy <- make_example_data(n_samples = 8)
    single_matrix <- dummy$omics$mRNA
    single_row_data <- dummy$row_data$mRNA

    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = list("Gene Expression" = single_matrix),
        row_data = list(single_row_data)
    )

    # does it work with just one experiment
    expect_s4_class(mae, "MultiAssayExperiment")
    expect_equal(length(experiments(mae)), 1)
    expect_true("Gene Expression" %in% names(experiments(mae)))
})

test_that("create_exposomicset generates row_data if not provided", {
    dummy <- make_example_data(n_samples = 6)

    mae <- create_exposomicset(
        codebook = dummy$codebook,
        exposure = dummy$exposure,
        omics = dummy$omics
        # no row_data provided
    )

    # ensure that some DataFrame is created when none is provided
    expect_s4_class(mae, "MultiAssayExperiment")
    expect_true(all(names(experiments(mae)) %in% c("mRNA", "proteomics")))
    expect_true(is(rowData(mae[["mRNA"]]), "DataFrame"))
})

test_that("create_exposomicset throws error on invalid input types", {
    dummy <- make_example_data()

    # error out if exposure is not a df
    expect_error(
        create_exposomicset(codebook = dummy$codebook, exposure = as.matrix(dummy$exposure), omics = dummy$omics),
        "must be a data frame"
    )

    # error out if omics is not a list or matrix
    expect_error(
        create_exposomicset(codebook = dummy$codebook, exposure = dummy$exposure, omics = 42),
        "must be a list or a single matrix"
    )

    # error out if the row data is not a list
    expect_error(
        create_exposomicset(codebook = dummy$codebook, exposure = dummy$exposure, omics = dummy$omics, row_data = "not_a_list"),
        "must be a list"
    )

    # error out if the rowdata does not match omics length
    expect_error(
        create_exposomicset(codebook = dummy$codebook, exposure = dummy$exposure, omics = dummy$omics, row_data = list(S4Vectors::DataFrame())),
        "Length of 'row_data' must match 'omics'"
    )
})
