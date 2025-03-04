filter_non_normal <- function(expomicset, p_thresh = 0.05) {
  require(MultiAssayExperiment)
  
  # Extract non-normal variables based on the p-value threshold
  non_normal_vars <- MultiAssayExperiment::metadata(expomicset)$transformation$norm_df$exposure[
    MultiAssayExperiment::metadata(expomicset)$transformation$norm_df$p.value < p_thresh
  ]
  
  message("Filtering out ", length(non_normal_vars), " non-normal exposure variables.")
  
  # Filter colData in the main object
  MultiAssayExperiment::colData(expomicset) <- MultiAssayExperiment::colData(expomicset)[
    , !colnames(MultiAssayExperiment::colData(expomicset)) %in% non_normal_vars, 
    drop = FALSE]
  
  # Filter colData in each experiment
  for (omics_name in names(experiments(expomicset))) {
    current_colData <- MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(expomicset)[[omics_name]])
    current_colData <- current_colData[
      , !colnames(current_colData) %in% non_normal_vars, 
      drop = FALSE]
    MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(expomicset)[[omics_name]]) <- current_colData
  }
  
  return(expomicset)
}
