#' Extract Feature Metadata from a MultiAssayExperiment
#'
#' Extracts feature-level metadata across all assays in a `MultiAssayExperiment` 
#' and returns a combined tibble.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#'
#' @details
#' This function:
#' - Iterates over all assays in the `MultiAssayExperiment`.
#' - Updates each assay’s sample metadata (`colData`) using `.update_assay_colData()`.
#' - Extracts feature-level metadata using `tidybulk::pivot_transcript()`.
#' - Combines results across assays into a single tibble, adding a `.exp_name` column.
#'
#' @return A tibble with feature metadata from all assays, with an added `.exp_name` column.
#'
#' @examples
#' \dontrun{
#' feature_data <- pivot_feature(expom)
#' }
#'
#' @export
pivot_feature <- function(expomicset){
  res <- lapply(names(MultiAssayExperiment::experiments(expomicset)),function(exp_name){
    exp <- .update_assay_colData(expomicset,exp_name) |>
      tidybulk::pivot_transcript()
  })
  
  names(res) <- names(MultiAssayExperiment::experiments(expomicset))
  
  res <- res |>
    dplyr::bind_rows(.id=".exp_name")
  return(res)
}