#' Perform Differential Abundance Analysis
#'
#' Runs differential abundance testing across all assays in a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param formula A model formula for differential testing.
#' @param abundance_col A character string specifying the assay data column to use. Default is `"counts"`.
#' @param minimum_counts An integer specifying the minimum count threshold for feature filtering. Default is `10`.
#' @param minimum_proportion A numeric specifying the minimum proportion of samples with nonzero counts. Default is `0.3`.
#' @param method A character string specifying the differential abundance method (`"limma_voom"`, `"DESeq2"`, etc.). Default is `"limma_voom"`.
#' @param contrasts A character vector specifying contrasts for differential testing. Default is `NULL`.
#' @param scaling_method A character string specifying normalization/scaling. Default is `"none"`.
#' @param action A character string specifying `"add"` (store results in metadata) or `"get"` (return results). Default is `"add"`.
#'
#' @details
#' This function:
#' - Iterates over all assays in `expomicset`.
#' - Updates sample metadata in each assay.
#' - Runs `.run_se_differential_abundance()` for differential testing.
#' - Filters results based on `minimum_counts` and `minimum_proportion`.
#' - Combines results across assays.
#' - Stores results in `metadata(expomicset)$differential_abundance` when `action="add"`.
#'
#' @return If `action="add"`, returns the updated `expomicset`.  
#' If `action="get"`, returns a `data.frame` with differential abundance results, including:
#' \item{feature}{Feature identifier.}
#' \item{exp_name}{Assay name.}
#' \item{logFC}{Log-fold change.}
#' \item{adj.P.Val}{Adjusted p-value.}
#'
#' @examples
#' \dontrun{
#' expom <- run_differential_abundance(
#'   expomicset = expom,
#'   formula = ~ condition,
#'   abundance_col = "counts",
#'   method = "limma_voom"
#' )
#' }
#'
#' @export
run_differential_abundance <- function(
    expomicset,
    formula,
    abundance_col = "counts",
    minimum_counts = 10,
    minimum_proportion = 0.3,
    method = "limma_voom",
    contrasts = NULL,
    scaling_method = "none",
    action="add"
) {
  
  message("Running differential abundance testing...")
  
  # Initialize a data frame to store results
  da_results_df <- list()
  
  # Iterate through assays in expomicset
  for (exp_name in names(MultiAssayExperiment::experiments(expomicset))) {
    message("Processing assay: ", exp_name)
    
    # Update assay with colData
    exp <- .update_assay_colData(expomicset, exp_name)
    
    # Run differential analysis using `.run_se_differential_abundance`
    res <- .run_se_differential_abundance(
      se = exp,
      formula = formula,
      abundance_col = abundance_col,
      method = method,
      scaling_method = scaling_method,
      min_counts = minimum_counts,
      min_proportion = minimum_proportion,
      contrasts = contrasts
    )
    
    # If results exist, append assay name and store them
    if (!is.null(res) && nrow(res) > 0) {
      res <- res |> 
        dplyr::mutate(exp_name = exp_name)
      da_results_df[[exp_name]] <- res
    } else {
      warning("No significant results found for assay: ", exp_name)
    }
  }
  
  # Combine results across assays
  final_results <- da_results_df |> 
    dplyr::bind_rows()
  
  message("Differential abundance testing completed.")
  
  if(action == "add") {
    MultiAssayExperiment::metadata(expomicset)$differential_abundance <- final_results
    return(expomicset)
  } else if (action == "get") {
    return(final_results)
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  return(expomicset)
}
