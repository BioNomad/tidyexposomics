correlate_factors_with_exposures <- function(
    expomicset, 
    exposures,
    action="add") {
  require(Hmisc)
  require(tidyverse)
  require(reshape2)
  
  message("Extracting factors from integration results...")
  
  # Get integration results
  integration_results <- MultiAssayExperiment::metadata(expomicset)$integration_results
  method_used <- integration_results$method
  
  # Extract factors based on integration method
  if (method_used == "MOFA") {
    message("Detected MOFA+ results, extracting factors...")
    factors <- MOFA2::get_factors(integration_results$result)
  } else if (method_used == "MCIA") {
    message("Detected MCIA results, extracting global scores...")
    factors <- nipalsMCIA::integration_results$result@global_scores
  } else {
    stop("Unsupported integration method: ", method_used)
  }
  
  if (is.null(factors)) {
    stop("No factors found in integration results.")
  }
  
  # Extract exposures from colData
  exp_data <- MultiAssayExperiment::colData(expomicset)[, exposures, drop = FALSE]
  
  # Ensure numeric data only
  factors <- as.data.frame(factors) |> 
    dplyr::mutate_all(as.numeric)
  exp_data <- as.data.frame(exp_data) |> 
    dplyr::mutate_all(as.numeric)
  
  # Ensure common samples
  common_samples <- intersect(rownames(factors), rownames(exp_data))
  if (length(common_samples) < 2) {
    stop("Not enough overlapping samples between factors and exposures for correlation.")
  }
  
  # Subset to common samples
  factors <- factors[common_samples, , drop = FALSE]
  exp_data <- exp_data[common_samples, , drop = FALSE]
  
  # Combine factor scores and exposures into a single matrix
  combined_df <- dplyr::bind_cols(factors, exp_data) %>%
    as.matrix()
  
  # Compute Spearman correlation using `rcorr()`
  message("Running Spearman correlation using rcorr()...")
  meta_cor <- Hmisc::rcorr(combined_df, type = "spearman")
  
  # Convert correlation and p-value matrices into long format
  meta_cor_df <- dplyr::inner_join(
    meta_cor$r |>
      as.data.frame() |>
      dplyr::rownames_to_column("var1") |>
      reshape2::melt(id.vars = "var1", 
           variable.name = "var2", 
           value.name = "correlation"),
    meta_cor$P |>
      as.data.frame() |>
      rownames_to_column("var1") |>
      reshape2::melt(id.vars = "var1", 
           variable.name = "var2", 
           value.name = "p.value"),
    by = c("var1", "var2")
  ) %>%
    # Keep only factor-exposure correlations (No exposure-exposure or factor-factor correlations)
    dplyr::filter(var1 %in% colnames(factors) &
                    var2 %in% colnames(exp_data)) |> 
    # Apply FDR correction
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) 
  
  message("Correlation analysis complete. Storing results...")
  
  if(action=="add"){
    # Store results in metadata
    MultiAssayExperiment::metadata(expomicset)$factor_exposure_correlations <- meta_cor_df
    return(expomicset)
  }else if (action=="get"){
    return(meta_cor_df)
  }else{
    
  }
}
