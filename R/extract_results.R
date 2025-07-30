#' Extract Results from `MultiAssayExperiment` Metadata
#'
#' Retrieves a specific analysis result
#' from the metadata slot of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param result A character string indicating which result to
#' extract from metadata. Must be one of:
#' `"codebook"`, `"quality_control"`, `"correlation"`, `"association"`,
#' `"differential_analysis"`, `"multiomics_integration"`, `"network"`,
#'  or `"enrichment"`.
#'
#' @return The corresponding result object stored in `metadata(expomicset)`,
#'  or `NULL` if not present.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 10,
#'   return_mae=TRUE
#'   )
#'
#' # extract results
#' res <- extract_results(
#'     expomicset=mae,
#'     result="codebook"
#' )
#'
#' @export
extract_results <- function(
    expomicset,
    result=c("codebook",
             "quality_control",
             "correlation",
             "association",
             "differential_analysis",
             "multiomics_integration",
             "network",
             "enrichment")){
  result <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck(result)

  return(result)

}
