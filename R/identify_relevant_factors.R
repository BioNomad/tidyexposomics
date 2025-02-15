identify_relevant_factors <- function(expOmicSet, outcome_var, categorical = FALSE, p_thresh = 0.05) {
  require(Hmisc)
  require(stats)
  require(broom)
  require(dplyr)
  require(purrr)
  require(tibble)
  
  message("Extracting factors from integration results...")
  
  # Get integration results
  integration_results <- expOmicSet@metadata$integration_results
  method_used <- integration_results$method
  
  # Extract factors based on integration method
  if (method_used == "MOFA") {
    message("Detected MOFA+ results, extracting factors...")
    factors <- get_factors(integration_results$result)
  } else if (method_used == "MCIA") {
    message("Detected MCIA results, extracting global scores...")
    factors <- integration_results$result@global_scores
  } else {
    stop("Unsupported integration method: ", method_used)
  }
  
  if (is.null(factors)) {
    stop("No factors found in integration results.")
  }
  
  # Extract outcome variable from colData
  outcome_data <- colData(expOmicSet) |> 
    as.data.frame()
  outcome_data <- outcome_data[, outcome_var, drop = FALSE]
  
  # Ensure common samples
  common_samples <- intersect(rownames(factors), rownames(outcome_data))
  if (length(common_samples) < 2) {
    stop("Not enough overlapping samples between factors and outcome for correlation.")
  }
  
  # Subset to common samples
  factors <- factors |> 
    as.data.frame() |> 
    (\(df) df |> filter(rownames(df) %in% common_samples))()
  outcome_data <- outcome_data |> 
    as.data.frame() |> 
    (\(df) df |> filter(rownames(df) %in% common_samples))()
  
  # Merge factors and outcome into a single data frame
  merged <- factors |> 
    as.data.frame() |> 
    rownames_to_column("id") |> 
    inner_join(outcome_data |> 
                 as.data.frame() |> 
                 rownames_to_column("id"),
               by = "id")
  
  message("Running correlation or Kruskal-Wallis test...")
  
  # Compute correlation/Kruskal-Wallis test using map2()
  results <- map2(colnames(factors), outcome_var, \(factor, outcome) {
    if (categorical) {
      kruskal.test(merged[[factor]] ~ merged[[outcome]]) |> tidy() |> 
        mutate(factor = factor, outcome = outcome)
    } else {
      cor.test(merged[[factor]], merged[[outcome]], method = "spearman", exact = FALSE) |> 
        tidy() |> mutate(factor = factor, outcome = outcome)
    }
  }) |> bind_rows()
  
  # Adjust for multiple testing (Benjamini-Hochberg FDR correction)
  #results <- results |> mutate(FDR = p.adjust(p.value, method = "fdr"))
  
  # Keep only significant results
  results <- results |> filter(p.value < p_thresh)
  
  message("Factors with significant association with outcome:")
  for (i in 1:nrow(results)) {
    message(paste0("Factor: ", results$factor[i], ", p-value: ", round(results$p.value[i],digits = 3)))
  }
  
  # Store significant factors in metadata
  expOmicSet@metadata$significant_factors <- results
  
  return(expOmicSet)
}
