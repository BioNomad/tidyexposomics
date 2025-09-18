#' Extract Feature Metadata from a MultiAssayExperiment
#'
#' Extracts feature-level metadata across all assays in a
#' `MultiAssayExperiment` and returns a combined tibble.
#'
#' @param exposomicset A `MultiAssayExperiment` object.
#'
#' @details
#' This function:
#' - Iterates over all assays in the `MultiAssayExperiment`.
#' - Updates each assay's sample metadata (`colData`) using
#' `.update_assay_colData()`.
#' - Extracts feature-level metadata using `pivot_transcript()` from `tidybulk`.
#' - Combines results across assays into a single tibble,
#' adding a `.exp_name` column.
#'
#' @return A tibble with feature metadata from all assays,
#' with an added `.exp_name` column.
#'
#'
#' @examples
#'
#' #' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # pivot experiment
#' feature_data <- mae |>
#'     pivot_feature()
#' @importMethodsFrom tidybulk pivot_transcript
#'
#' @export
pivot_feature <- function(exposomicset) {
    res <- lapply(
        names(MultiAssayExperiment::experiments(exposomicset)), function(exp_name) {
            exp <- .update_assay_colData(exposomicset, exp_name) |>
                pivot_transcript()
        }
    )

    names(res) <- names(MultiAssayExperiment::experiments(exposomicset))

    res <- res |>
        dplyr::bind_rows(.id = ".exp_name")
    return(res)
}
