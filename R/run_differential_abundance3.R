run_differential_abundance <- function(
    expOmicSet,
    formula,
    abundance_col = "counts",
    minimum_counts = 10,
    minimum_proportion = 0.3,
    method = "limma_voom",
    contrasts = NULL,
    scaling_method = "none"
) {
  require(tidybulk)
  require(MultiAssayExperiment)
  require(tidyverse)
  
  message("Running differential abundance testing...")
  
  # Initialize a data frame to store results
  da_results_df <- list()
  
  # Iterate through assays in expOmicSet
  for (exp_name in names(experiments(expOmicSet))) {
    message("Processing assay: ", exp_name)
    
    # Update assay with colData
    assay <- .update_assay_colData(expOmicSet, exp_name)
    
    # Run differential analysis using `.run_se_differential_abundance`
    res <- .run_se_differential_abundance(
      exp = assay,
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
  
  # Save results in metadata if available
  if (nrow(final_results) > 0) {
    metadata(expOmicSet)$differential_abundance <- final_results |> drop_na()
  } else {
    warning("No differential abundance results found across assays.")
  }
  
  message("Differential abundance testing completed.")
  return(expOmicSet)
}
