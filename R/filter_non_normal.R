filter_non_normal <- function(expOmicSet, p_thresh = 0.05) {
  # Extract non-normal variables based on the p-value threshold
  non_normal_vars <- metadata(expOmicSet)$transformation$norm_df$exposure[
    metadata(expOmicSet)$transformation$norm_df$p.value < p_thresh
  ]
  
  message("Filtering out ", length(non_normal_vars), " non-normal exposure variables.")
  
  # Filter colData in the main object
  colData(expOmicSet) <- colData(expOmicSet)[, !colnames(colData(expOmicSet)) %in% non_normal_vars, drop = FALSE]
  
  # Filter colData in each experiment
  for (omics_name in names(experiments(expOmicSet))) {
    current_colData <- colData(experiments(expOmicSet)[[omics_name]])
    current_colData <- current_colData[, !colnames(current_colData) %in% non_normal_vars, drop = FALSE]
    colData(experiments(expOmicSet)[[omics_name]]) <- current_colData
  }
  
  return(expOmicSet)
}
