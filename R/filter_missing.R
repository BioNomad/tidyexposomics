# filter_missing <- function(expomicset, 
#                            na_thresh = 20, 
#                            na_plot_thresh = 5) {
#   require(naniar)
#   require(tidyverse)
#   
#   # Helper function to identify variables/features to exclude
#   get_vars_to_exclude <- function(data, thresh) {
#     miss_summary <- naniar::miss_var_summary(data)
#     list(
#       to_exclude = miss_summary %>% filter(pct_miss > thresh) %>% pull(variable),
#       summary = miss_summary
#     )
#   }
#   
#   # Handle metadata (colData)
#   exposure <- as.data.frame(colData(expomicset))
#   col_exclude <- get_vars_to_exclude(exposure, na_thresh)
#   colData(expomicset) <- as(colData(expomicset)[
#     , !colnames(colData(expomicset)) %in% 
#       col_exclude$to_exclude, drop = FALSE], 
#     "DataFrame")
#   
#   # QC plot for metadata
#   metadata(expomicset)$na_qc <- list(
#     exposure = list(
#       vars_to_exclude_exposure_sum = col_exclude$summary %>% filter(pct_miss > na_thresh),
#       na_exposure_qc_plot = naniar::gg_miss_var(
#         exposure %>% dplyr::select(all_of(col_exclude$to_exclude))
#       )
#     )
#   )
#   
#   message("Missing Data Filter threshold: ", na_thresh, "%")
#   message("Filtered metadata variables: ", paste(col_exclude$to_exclude, collapse = ", "))
#   
#   # Handle omics experiments
#   for (omics_name in names(experiments(expomicset))) {
#     experiment <- experiments(expomicset)[[omics_name]]
#     original_assay <- assays(experiment)[[1]]
#     
#     # Transpose and convert assay to data.frame for processing
#     assay_data <- t(original_assay) %>% as.data.frame()
#     
#     # Identify features (rows) to keep based on missingness
#     row_exclude <- get_vars_to_exclude(assay_data, na_thresh)
#     
#     # Filter colData variables with high missingness in this experiment
#     filtered_colData_exp <- as(colData(experiment)[
#       , !colnames(colData(experiment)) %in% 
#         col_exclude$to_exclude, drop = FALSE], 
#       "DataFrame")
#     
#     # Update the experiment
#     colData(experiment) <- filtered_colData_exp
#     experiment <- experiment[!rownames(experiment) %in% row_exclude$to_exclude,]
#     
#     experiments(expomicset)[[omics_name]] <- experiment
#     
#     # QC plot for omics data
#     metadata(expomicset)$na_qc[[omics_name]] <- list(
#       vars_to_exclude_omics_sum = row_exclude$summary %>% filter(pct_miss > na_thresh),
#       na_omics_qc_plot = naniar::gg_miss_var(assay_data %>% dplyr::select(all_of(row_exclude$to_exclude)))
#     )
#     
#     message("Filtered rows with high missingness in ", omics_name, ": ", length(row_exclude$to_exclude))
#   }
#   
#   return(expomicset)
# }


filter_missing <- function(expomicset, 
                           na_thresh = 20, 
                           na_plot_thresh = 5) {
  require(naniar)
  require(tidyverse)
  
  # Helper function to identify variables/features to exclude
  get_vars_to_exclude <- function(data, thresh) {
    miss_summary <- naniar::miss_var_summary(data)
    list(
      to_exclude = miss_summary %>% filter(pct_miss > thresh) %>% pull(variable),
      summary = miss_summary
    )
  }
  
  # Handle metadata (colData)
  exposure <- as.data.frame(colData(expomicset))
  col_exclude <- get_vars_to_exclude(exposure, na_thresh)
  colData(expomicset) <- as(colData(expomicset)[
    , !colnames(colData(expomicset)) %in% 
      col_exclude$to_exclude, drop = FALSE], 
    "DataFrame")
  
  # QC plot for metadata
  metadata(expomicset)$na_qc <- list(
    exposure = list(
      vars_to_exclude_exposure_sum = col_exclude$summary %>% filter(pct_miss > na_thresh),
      all_var_sum = col_exclude$summary %>% filter(pct_miss > 0),
      na_exposure_qc_plot = naniar::gg_miss_var(
        exposure %>% dplyr::select(all_of(col_exclude$to_exclude))
      )
    )
  )
  
  message("Missing Data Filter threshold: ", na_thresh, "%")
  message("Filtered metadata variables: ", paste(col_exclude$to_exclude, collapse = ", "))
  
  # Handle omics experiments
  for (omics_name in names(experiments(expomicset))) {
    experiment <- experiments(expomicset)[[omics_name]]
    original_assay <- assays(experiment)[[1]]
    
    # Transpose and convert assay to data.frame for processing
    assay_data <- as.data.frame(t(original_assay))
    
    # Identify features (rows) to exclude based on missingness
    row_exclude <- get_vars_to_exclude(assay_data, na_thresh)
    
    # Filter rows with high missingness
    experiment <- experiment[!rownames(experiment) %in% row_exclude$to_exclude,]

    experiments(expomicset)[[omics_name]] <- experiment
    
    # QC plot for omics data
    metadata(expomicset)$na_qc[[omics_name]] <- list(
      vars_to_exclude_omics_sum = row_exclude$summary %>% filter(pct_miss > na_thresh),
      all_var_sum = row_exclude$summary %>% filter(pct_miss > 0),
      na_omics_qc_plot = naniar::gg_miss_var(assay_data %>% dplyr::select(all_of(row_exclude$to_exclude)))
    )
    
    message("Filtered rows with high missingness in ", omics_name, ": ", length(row_exclude$to_exclude))
  }
  
  return(expomicset)
}

