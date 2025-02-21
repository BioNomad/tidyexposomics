extract_top_factor_features <- function(
    expOmicSet, 
    factors, 
    method = "percentile", 
    percentile = 0.9, 
    threshold = 0.3) {
  require(dplyr)
  require(tidyr)
  require(purrr)
  
  message("Extracting top contributing features for specified factors...")
  
  # Get integration results
  integration_results <- metadata(expOmicSet)$integration_results
  method_used <- integration_results$method
  
  # Extract factor loadings based on method
  if (method_used == "MOFA") {
    message("Using MOFA+ factor loadings...")
    loadings <- get_weights(integration_results$result)  
  } else if (method_used == "MCIA") {
    message("Using MCIA block loadings...")
    loadings <- integration_results$result@block_loadings  
  } else {
    stop("Unsupported integration method: ", method_used)
  }
  
  # Convert factor loadings to long format
  loadings_df <- loadings |> 
    map2(names(loadings), \(df, exp_name) {
      df <- as.data.frame(df)  
      df$exp_name <- exp_name  
      df$feature <- rownames(df)  
      rownames(df) <- NULL  
      return(df) 
    }) |> bind_rows() |> 
    pivot_longer(
      -c(feature, exp_name),
      names_to = "factor",
      values_to = "loading")
  
  # Ensure factor names are treated as characters and match MCIA format
  factors <- as.character(factors)
  loadings_df <- loadings_df |>
    mutate(factor = as.character(factor))
  
  # Filter for user-specified factors
  loadings_df <- loadings_df |> 
    filter(factor %in% factors)
  
  # Apply thresholding method
  if (method == "percentile") {
    message(
      "Applying percentile-based filtering (",
      percentile * 100,
      "%)...")
    
    factor_thresholds <- loadings_df |> 
      group_by(factor) |> 
      dplyr::summarize(
        threshold = quantile(abs(loading),
                             percentile,
                             na.rm = TRUE),
                .groups = "drop")
  } else if (method == "threshold") {
    message(
      "Applying raw threshold-based filtering (>|",
      threshold,
      "|)...")
    factor_thresholds <- tibble(
      factor = unique(loadings_df$factor),
      threshold = threshold)
  } else {
    stop("Invalid method. Choose 'percentile' or 'threshold'.")
  }
  
  # Merge computed thresholds with loadings
  loadings_df <- left_join(
    loadings_df, 
    factor_thresholds, 
    by = "factor")
  
  # Apply filtering
  filtered_features <- loadings_df |> 
    filter(abs(loading) > threshold) |> 
    dplyr::select(feature,
                  factor,
                  loading,
                  exp_name) |>
    distinct()
  
  # Store selected features
  metadata(expOmicSet)$top_factor_features <- filtered_features
  
  message("Selected ",
          nrow(filtered_features), 
          " features contributing to specified factors.")
  
  return(expOmicSet)
}
