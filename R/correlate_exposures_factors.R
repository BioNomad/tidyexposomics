#' Correlate Factor Features with Exposures
#'
#' Performs correlation analysis between factor-contributing features and exposure variables.
#' Supports batch processing for large datasets and multiple testing correction.
#'
#' @param expomicset A `MultiAssayExperiment` object containing factor feature information.
#' @param exposure_cols An optional character vector of exposure variable names. If `NULL`, all numeric exposures are used. Default is `NULL`.
#' @param correlation_method A character string specifying the correlation method (`"spearman"`, `"pearson"`, or `"kendall"`). Default is `"spearman"`.
#' @param correlation_cutoff A numeric threshold for absolute correlation values to retain results. Default is `0.3`.
#' @param cor_pval_column A character string specifying the column name for correlation p-values. Default is `"p.value"`.
#' @param pval_cutoff A numeric threshold for statistical significance in correlation analysis. Default is `0.05`.
#' @param batch_size An integer specifying the number of features to process in each batch. Default is `1500`.
#' @param action A character string indicating whether to return results (`"get"`) or add them to metadata (`"add"`). Default is `"add"`.
#'
#' @details
#' The function extracts factor-contributing features from `metadata(expomicset)`, filters them based on experiment type, 
#' and correlates them with numeric exposure variables in `colData(expomicset)`. Correlations are computed in batches 
#' using the specified method, and p-values are adjusted using FDR correction.
#'
#' @return If `action = "add"`, returns the modified `expomicset` with correlation results stored in metadata.
#' If `action = "get"`, returns a data frame containing:
#' \item{exposure}{The exposure variable tested.}
#' \item{feature}{The factor feature tested.}
#' \item{estimate}{The correlation coefficient.}
#' \item{p.value}{The correlation p-value.}
#' \item{FDR}{The adjusted p-value (false discovery rate).}
#' \item{exp_name}{The experiment from which the feature originated.}
#'
#' @examples
#' \dontrun{
#' results <- correlate_exposures_factors(
#'   expomicset = expom,
#'   exposure_cols = c("PM2.5", "NO2"),
#'   correlation_method = "spearman",
#'   action = "get"
#' )
#' }
#'
#' @export
correlate_exposures_factors <- function(
    expomicset,
    exposure_cols = NULL,  
    correlation_method = "spearman",
    correlation_cutoff = 0.3,        
    cor_pval_column = "p.value",      
    pval_cutoff = 0.05,
    batch_size = 1500,
    action = "add"
) {

  message("Starting correlation analysis between factor features and exposures...")
  
  # Extract factor-contributing features
  factor_features <- MultiAssayExperiment::metadata(expomicset)$top_factor_features
  if (is.null(factor_features)) {
    stop("No factor features found in metadata.")
  }
  
  # Get numeric exposure variables
  numeric_exposures <- colnames(MultiAssayExperiment::colData(expomicset))
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
      dplyr::filter(exp_name == experiment_name) |> 
      dplyr::pull(feature) |> 
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
      correlation_results[[experiment_name]] <- batch_results |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(exp_name = experiment_name) 
    }
  }
  
  combined_results <- correlation_results |> 
    dplyr::bind_rows() |> 
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |> 
    dplyr::left_join(
      MultiAssayExperiment::metadata(expomicset)$var_info,
      by=c("exposure"="variable"))
  
  if (nrow(combined_results) == 0) {
    warning("No significant correlations found in any experiment.")
    return(expomicset)
  }
  
  if(action == "add") {
    # Save to metadata
    MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation <- combined_results
    return(expomicset)
  } else if (action == "get") {
    return(combined_results)
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  message("Factor feature-exposure correlation analysis completed.")
  
}
