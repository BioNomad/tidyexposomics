correlate_exposures_with_factors <- function(
    expomicset,
    exposure_cols = NULL,  
    correlation_method = "spearman",
    correlation_cutoff = 0.3,        
    cor_pval_column = "p.value",      
    pval_cutoff = 0.05,
    batch_size = 1500,
    action = "add"
) {
  require(tidyverse)
  require(MultiAssayExperiment)
  
  message("Starting correlation analysis between factor features and exposures...")
  
  # Extract factor-contributing features
  factor_features <- metadata(expomicset)$top_factor_features
  if (is.null(factor_features)) {
    stop("No factor features found in metadata.")
  }
  
  # Get numeric exposure variables
  numeric_exposures <- colnames(colData(expomicset))
  if (!is.null(exposure_cols)) {
    numeric_exposures <- intersect(numeric_exposures, exposure_cols)
  }
  if (length(numeric_exposures) == 0) {
    stop("No numeric exposure variables found in colData.")
  }
  
  correlation_results <- list()
  
  for (experiment_name in unique(factor_features$exp_name)) {
    message("Processing experiment: ", experiment_name)
    
    # Extract relevant factor features for this assay
    selected_features <- factor_features |> 
      filter(exp_name == experiment_name) |> 
      pull(feature) |> 
      unique()
    
    if (length(selected_features) == 0) {
      warning("No relevant factor features found for ", 
              experiment_name,
              ", skipping.")
      next
    }
    
    # Extract and update SummarizedExperiment
    se <- .update_assay_colData(expomicset, experiment_name)
    
    # **Filter SummarizedExperiment to Selected Features**
    se <- se[rownames(se) %in% selected_features, , drop = FALSE]
    
    if (nrow(se) == 0) {
      warning(
        "No valid factor features left in assay after filtering for ",
        experiment_name,
        ", skipping.")
      next
    }
    
    # **Batch Processing**
    feature_batches <- split(
      selected_features, 
      ceiling(seq_along(selected_features) / batch_size))
    
    batch_results <- list()
    
    batch_index <- 1
    
    for (batch in feature_batches) {
      message("  - Processing batch ",
              batch_index, 
              " of ",
              length(feature_batches), 
              " (", length(batch),
              " features)...")
      
      batch_index <- batch_index + 1
      
      # **Subset Data to Only Batch Features**
      se_batch <- se[rownames(se) %in% batch, , drop = FALSE]
      
      # Perform correlation analysis
      batch_result <- .correlate_se_with_coldata(
        se = se_batch,
        exposure_cols = numeric_exposures,
        correlation_method = correlation_method,
        correlation_cutoff = correlation_cutoff,
        cor_pval_column = cor_pval_column,
        pval_cutoff = pval_cutoff
      )
      
      if (nrow(batch_result) > 0) {
        batch_results[[length(batch_results) + 1]] <- batch_result
      }
    }
    
    if (length(batch_results) > 0) {
      correlation_results[[experiment_name]] <- bind_rows(batch_results) |> 
        mutate(exp_name = experiment_name) 
    }
  }
  
  combined_results <- bind_rows(correlation_results) |> 
    mutate(FDR = p.adjust(p.value, method = "fdr")) |> 
    left_join(expomicset@metadata$var_info,
              by=c("exposure"="variable"))
  
  if (nrow(combined_results) == 0) {
    warning("No significant correlations found in any experiment.")
    return(expomicset)
  }
  
  if(action == "add") {
    # Save to metadata
    metadata(expomicset)$omics_exposure_factor_correlation <- combined_results
    return(expomicset)
  } else if (action == "get") {
    return(combined_results)
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  message("Factor feature-exposure correlation analysis completed.")
  
}
