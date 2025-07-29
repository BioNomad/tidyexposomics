#' Extract Results from `MultiAssayExperiment` Metadata
#'
#' Retrieves a specific analysis result (e.g., quality control, correlation, association, etc.)
#' from the metadata slot of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param result A character string indicating which result to extract from metadata. Must be one of:
#' `"codebook"`, `"quality_control"`, `"correlation"`, `"association"`,
#' `"differential_analysis"`, `"multiomics_integration"`, `"network"`, or `"enrichment"`.
#'
#' @return The corresponding result object stored in `metadata(expomicset)`, or `NULL` if not present.
#'
#' @examples
#' \dontrun{
#' extract_results(expomicset, result = "association")
#' extract_results(expomicset, result = "multiomics_integration")
#' }
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
