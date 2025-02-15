exposure_omic_association <- function(
    expOmicSet, 
    exposures, 
    da_p_value_column = "adj.P.Val",
    da_logFC_column = "logFC",
    da_threshold = 0.05,
    da_logFC_threshold = log2(1.5),
    confounders = NULL, 
    family = "gaussian", 
    correction_method = "fdr"
) {
  library(tidyverse)
  library(broom)
  
  # Extract differential abundance results
  message("Extracting differential abundance results...")
  da_res <- metadata(expOmicSet)$differential_abundance |> 
    filter(!!sym(da_p_value_column) < da_threshold) |> 
    filter(abs(!!sym(da_logFC_column)) > da_logFC_threshold)
  
  # Validate differential abundance results
  if (is.null(da_res) || nrow(da_res) == 0) {
    stop("No differential abundance results found in metadata.")
  }
  
  # Define variables with non-zero variance
  non_zero <- colData(expOmicSet) |> 
    as.data.frame() |> 
    select(where(is.numeric)) |> 
    map_dbl(var) |> 
    keep(~ . > 0) |> 
    names()
  
  # Initialize results list
  results <- list()
  
  # Iterate over experiments
  message("Performing exposure-omic association for each experiment...")
  for (experiment_name in unique(da_res$assay_name)) {
    message(" - Analyzing experiment: ", experiment_name)
    
    # Update the assay with colData for the current experiment
    experiment <- .update_assay_colData(expOmicSet, experiment_name)
    
    # Filter DEGs for the current experiment
    da_features <- da_res |> 
      filter(assay_name == experiment_name) |> 
      pull(molecular_feature)
    
    # Generate a random merge column name to ensure it does not overlap with existing column names
    merge_col <- stringi::stri_rand_strings(1, 5)
    
    # Process colData and merge with assay data
    merged_data <- colData(experiment) |> 
      as.data.frame() |> 
      mutate(across(all_of(non_zero), ~ as.numeric(scale(.)))) |>  # Scale numeric variables
      rownames_to_column(var = merge_col) |> 
      inner_join(
        assays(experiment)[[1]] |> 
          t() |> 
          as.data.frame() |> 
          select(all_of(da_features)) |> 
          mutate(across(everything(), ~ as.numeric(scale(log2(. + 1))))) |> 
          rownames_to_column(var = merge_col),
        by = merge_col
      )
    
    # Define exposure columns
    exposure_cols <- colnames(merged_data) |> 
      intersect(colnames(colData(expOmicSet)))
    
    exposure_cols <- exposure_cols |> 
      intersect(exposures)
    
    # For each DEG, test association with exposures
    exp_res <- map_dfr(da_features, function(feature) {
      map_dfr(exposure_cols, function(exposure) {
        # Construct formula
        formula <- as.formula(
          paste(feature, "~", exposure, if (!is.null(confounders)) paste("+", paste(confounders, collapse = " +")) else "")
        )
        
        # Fit the model
        tryCatch(
          glm(formula, data = merged_data, family = family) |> 
            tidy() |> 
            filter(term == exposure) |>  # Focus on the exposure term
            mutate(feature = feature, experiment = experiment_name),
          error = function(e) {
            message("Model fitting failed for feature ", feature, " and exposure ", exposure, ": ", e$message)
            return(NULL)
          }
        )
      })
    })
    
    # Append results for the current experiment
    results[[experiment_name]] <- exp_res
  }
  
  # Combine results across experiments
  message("Combining results...")
  results_df <- bind_rows(results)
  
  # Adjust p-values for multiple testing
  message("Applying multiple testing correction...")
  results_df <- results_df |> 
    mutate(p_adjusted = p.adjust(p.value, method = correction_method)) |> 
    dplyr::rename("exposure"=term)
  
  # Reorder columns for clarity
  results_df <- results_df |> 
    select(
      experiment, feature, exposure, estimate, std.error, statistic, p.value, p_adjusted, everything()
    ) |> 
    left_join(metadata(expOmicSet)$var_info, by = c("exposure" = "variable"))
  
  # Save results in metadata
  message("Saving results to expOmicSet metadata...")
  metadata(expOmicSet)$exposure_omic_association <- results_df
  
  message("Exposure-omic association analysis completed.")
  return(expOmicSet)
}
