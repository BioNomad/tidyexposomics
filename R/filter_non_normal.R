#' Filter Non-Normal Exposure Variables
#'
#' Removes exposure variables that deviate significantly from a normal distribution based on 
#' normality test results stored in metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure and omics data.
#' @param p_thresh A numeric value specifying the p-value threshold for normality. Variables with `p.value < p_thresh` are removed. Default is `0.05`.
#'
#' @details
#' The function identifies exposure variables that fail a normality test using `metadata(expomicset)$transformation$norm_df`. 
#' - Exposure variables with `p.value < p_thresh` are removed from `colData(expomicset)`.
#' - The same filtering is applied to `colData` in each experiment within `experiments(expomicset)`.
#'
#' @return A `MultiAssayExperiment` object with non-normal exposure variables removed.
#'
#' @examples
#' \dontrun{
#' filtered_expom <- filter_non_normal(
#'   expomicset = expom,
#'   p_thresh = 0.05)
#' }
#'
#' @export
filter_non_normal <- function(expomicset, 
                              p_thresh = 0.05) {
  
  # Extract non-normal variables based on the p-value threshold
  non_normal_vars <- MultiAssayExperiment::metadata(expomicset)$transformation$norm_df$exposure[
    MultiAssayExperiment::metadata(expomicset)$transformation$norm_df$p.value < p_thresh
  ]
  
  message("Filtering out ", length(non_normal_vars), " non-normal exposure variables.")
  
  # Filter colData in the main object
  MultiAssayExperiment::colData(expomicset) <- MultiAssayExperiment::colData(expomicset)[
    , !colnames(MultiAssayExperiment::colData(expomicset)) %in% non_normal_vars, 
    drop = FALSE]
  
  # Filter colData in each experiment
  for (omics_name in names(MultiAssayExperiment::experiments(expomicset))) {
    current_colData <- MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(expomicset)[[omics_name]])
    current_colData <- current_colData[
      , !colnames(current_colData) %in% non_normal_vars, 
      drop = FALSE]
    MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(expomicset)[[omics_name]]) <- current_colData
  }
  
  return(expomicset)
}
