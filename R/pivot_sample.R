#' Extract Sample Metadata from MultiAssayExperiment or SummarizedExperiment
#'
#' Extracts and formats sample-level metadata (`colData`) from a `MultiAssayExperiment`
#' or `SummarizedExperiment` object.
#'
#' @param x A `MultiAssayExperiment` or `SummarizedExperiment` object.
#' @param ... Additional arguments passed to `tidybulk::pivot_sample()` for `SummarizedExperiment` objects.
#'
#' @details
#' This function:
#' - Extracts **sample metadata** from `MultiAssayExperiment` using `colData()`, converting it to a tibble.
#' - Calls `tidybulk::pivot_sample()` when applied to a `SummarizedExperiment` object.
#' - **Error Handling**: Returns an error if `x` is not a `MultiAssayExperiment` or `SummarizedExperiment`.
#'
#' @return A tibble containing sample metadata with an added `.sample` column.
#'
#' @examples
#' \dontrun{
#' sample_data <- pivot_sample(expom)
#' }
#'
#' @export
pivot_sample <- function(x, ...) {
  if (inherits(x, "MultiAssayExperiment")) {
    # Call your custom function
    MultiAssayExperiment::colData(x) |>
      as.data.frame() |>
      tibble::rownames_to_column(".sample") |>
      dplyr::as_tibble()
  } else if (inherits(x, "SummarizedExperiment")) {
    # Call the pivot_sample() from the other package
    tidybulk::pivot_sample(x, ...)
  } else {
    stop("Error: pivot_sample() only supports MultiAssayExperiment and SummarizedExperiment objects.")
  }
}

