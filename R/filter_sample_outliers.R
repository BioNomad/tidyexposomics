#' Filter Sample Outliers
#'
#' Removes sample outliers from a `MultiAssayExperiment` object
#' based on PCA analysis.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing
#' omics and exposure data.
#' @param outliers An optional character vector specifying
#' sample names to be removed.
#' If `NULL`, the function uses outliers identified in
#' `metadata(exposomicset)$pca$outliers`. Default is `NULL`.
#'
#' @details
#' The function checks for the presence of PCA results in
#' `metadata(exposomicset)`. If `outliers` is not provided,
#' it retrieves precomputed outliers from `metadata(exposomicset)$pca$outliers`.
#' The identified samples are removed
#' from the dataset.
#'
#' @return A `MultiAssayExperiment` object with the specified outliers removed.
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run PCA
#' mae <- mae |>
#'     run_pca()
#'
#' # filter outliers if present
#' mae <- mae |>
#'     filter_sample_outliers()
#'
#' @export
filter_sample_outliers <- function(
    exposomicset,
    outliers = NULL) {
    # Check if the input is a MultiAssayExperiment object
    if (!"pca" %in% names(
        MultiAssayExperiment::metadata(exposomicset)$quality_control
    )) {
        stop("PCA not performed on the data. Please run 'pca_analysis' first.")
    }

    # Check if outliers are provided
    if (is.null(outliers)) {
        outliers <- exposomicset |>
            MultiAssayExperiment::metadata() |>
            purrr::pluck("quality_control") |>
            purrr::pluck("pca") |>
            purrr::pluck("outliers")
    } else {
        outliers <- outliers
    }

    message("Removing outliers: ", paste(outliers, collapse = ", "))

    # Subset by ColData to properly filter outliers
    exposomicset <- MultiAssayExperiment::subsetByColData(
        exposomicset,
        !(rownames(SummarizedExperiment::colData(exposomicset)) %in% outliers)
    )


    # Add analysis steps taken to metadata
    step_record <- list(filter_sample_outliers = list(
        timestamp = Sys.time(),
        params = list(),
        notes = paste("Outliers: ", paste(outliers, collapse = ", "), collapse = "")
    ))

    MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
        MultiAssayExperiment::metadata(exposomicset)$summary$steps,
        step_record
    )

    # exposomicset <- exposomicset[,!colnames(exposomicset) %in% outliers]
    return(exposomicset)
}
