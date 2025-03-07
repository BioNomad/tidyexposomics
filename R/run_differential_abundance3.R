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
  require(tidybulk)
  require(MultiAssayExperiment)
  require(tidyverse)
  
  message("Running differential abundance testing...")
  
  # Initialize a data frame to store results
  da_results_df <- list()
  
  # Iterate through assays in expomicset
  for (exp_name in names(experiments(expomicset))) {
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
      res <- res |> mutate(exp_name = exp_name)
      da_results_df[[exp_name]] <- res
    } else {
      warning("No significant results found for assay: ", exp_name)
    }
  }
  
  # Combine results across assays
  final_results <- bind_rows(da_results_df)
  
  message("Differential abundance testing completed.")
  
  if(action == "add") {
    metadata(expomicset)$differential_abundance <- final_results
    return(expomicset)
  } else if (action == "get") {
    return(final_results)
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  return(expomicset)
}
