perform_exwas <- function(expomicset, 
                          exposures, 
                          outcome,
                          covariates = NULL,
                          family = "gaussian",
                          correction_method = "fdr",
                          action="add") {
  library(tidyverse)
  library(broom)
  
  # Extract and preprocess colData
  message("Extracting exposure data...")
  data <- colData(expomicset)  |> 
    as.data.frame() |> 
    mutate_if(is.numeric, ~ scale(.))  
  
  # Initialize results list
  results <- list()
  
  # Iterate over each exposure
  message("Performing ExWAS...")
  for (exposure in exposures) {
    message(" - Analyzing exposure: ", exposure)
    
    # Construct formula
    formula <- as.formula(
      paste(outcome, "~", exposure, if (!is.null(covariates)) paste("+", paste(covariates, collapse = " +")) else "")
    )
    
    # Fit the generalized linear model
    model <- tryCatch(
      glm(formula, data = data, family = family),
      error = function(e) {
        message("   * Model fitting failed for ", exposure, ": ", e$message)
        return(NULL)
      }
    )
    
    # Skip if model failed to fit
    if (is.null(model)) next
    
    # Extract model statistics
    model_summary <- tidy(model) |>
      filter(term == exposure)  # Focus on the exposure term
    
    # Add exposure name
    model_summary <- model_summary |>
      mutate(exposure = exposure)
    
    # Append to results
    results[[exposure]] <- model_summary
  }
  
  # Combine results into a single data frame
  message("Combining results...")
  results_df <- bind_rows(results)
  
  # Adjust p-values for multiple testing
  message("Applying multiple testing correction...")
  results_df <- results_df |>
    mutate(p_adjusted = p.adjust(p.value, method = correction_method)) |> 
    mutate(depend = outcome)
  
  # Reorder columns for clarity
  results_df <- results_df |>
    dplyr::select(exposure, depend, estimate, std.error, statistic, p.value, p_adjusted, everything()) |> 
    dplyr::rename("outcome" = "depend") |> 
    left_join(metadata(expomicset)$var_info,
              by = c("exposure" = "variable"))
  
  message("ExWAS analysis completed.")
  
  if(action=="add"){
    # Save results in metadata
    metadata(expomicset)$exwas <- list(
      results_df=results_df,
      covariates=covariates)
    
    return(expomicset)
  }else if (action=="get"){
    return(list(
      results_df=results_df,
      covariates=covariates
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
