#' Correlate Differentially Expressed Genes (DEGs) with Exposures
#'
#' Performs correlation analysis between DEGs and exposure variables, with optional
#' filtering based on sensitivity analysis results. Supports batch processing for
#' large datasets and multiple testing correction.
#'
#' @param expomicset A `MultiAssayExperiment` object containing differential expression results.
#' @param exposure_cols An optional character vector of exposure variable names. If `NULL`, all numeric exposures are used. Default is `NULL`.
#' @param robust A logical value indicating whether to use sensitivity analysis to filter DEGs. Default is `TRUE`.
#' @param score_thresh An optional numeric threshold for feature stability scores in sensitivity analysis. If `NULL`, the default threshold from metadata is used. Default is `NULL`.
#' @param correlation_method A character string specifying the correlation method (`"spearman"`, `"pearson"`, or `"kendall"`). Default is `"spearman"`.
#' @param correlation_cutoff A numeric threshold for absolute correlation values to retain results. Default is `0.3`.
#' @param cor_pval_column A character string specifying the column name for correlation p-values. Default is `"p.value"`.
#' @param pval_cutoff A numeric threshold for statistical significance in correlation analysis. Default is `0.05`.
#' @param batch_size An integer specifying the number of features to process in each batch. Default is `1500`.
#' @param deg_pval_col A character string specifying the column name for DEG p-values. Default is `"adj.P.Val"`.
#' @param deg_logfc_col A character string specifying the column name for DEG log fold changes. Default is `"logFC"`.
#' @param deg_pval_thresh A numeric threshold for filtering DEGs based on p-values. Default is `0.05`.
#' @param deg_logfc_thresh A numeric threshold for filtering DEGs based on absolute log fold change. Default is `log2(1.5)`.
#' @param action A character string indicating whether to return results (`"get"`) or add them to metadata (`"add"`). Default is `"add"`.
#'
#' @details
#' The function extracts DEGs from `metadata(expomicset)`, filters them based on p-value and log fold change thresholds,
#' and correlates them with numeric exposure variables in `colData(expomicset)`. If `robust = TRUE`, only stable DEGs
#' from sensitivity analysis are used. Correlations are computed in batches using the specified method, and p-values are
#' adjusted using the Benjamini-Hochberg FDR correction.
#'
#' @return If `action = "add"`, returns the modified `expomicset` with correlation results stored in metadata.
#' If `action = "get"`, returns a data frame containing:
#' \item{exposure}{The exposure variable tested.}
#' \item{feature}{The DEG tested.}
#' \item{estimate}{The correlation coefficient.}
#' \item{p.value}{The correlation p-value.}
#' \item{FDR}{The adjusted p-value (false discovery rate).}
#' \item{exp_name}{The experiment from which the DEG originated.}
#'
#' @examples
#' \dontrun{
#' results <- correlate_exposures_degs(
#'   expomicset = expom,
#'   exposure_cols = c("PM2.5", "NO2"),
#'   robust = TRUE,
#'   correlation_method = "spearman",
#'   action = "get"
#' )
#' }
#'
#' @export
correlate_exposures_degs <- function(
    expomicset,
    exposure_cols = NULL,
    robust = TRUE,
    score_thresh = NULL,
    correlation_method = "spearman",
    correlation_cutoff = 0.3,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05,
    batch_size = 1500,
    deg_pval_col = "adj.P.Val",
    deg_logfc_col = "logFC",
    deg_pval_thresh = 0.05,
    deg_logfc_thresh = log2(1.5),
    action="add"
) {
  message("Starting correlation analysis between DEGs and exposures...")

  # Extract and filter DEGs
  da_results <- MultiAssayExperiment::metadata(expomicset)$differential_abundance
  if (is.null(da_results)) {
    stop("No differential abundance results found in metadata.")
  }

  # Stop if robust is true but sensitivity analysis was not run
  if(robust){
    if(!"sensitivity_analysis" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `run_sensitivity_analysis()` first.")
    }
  }

  # Filter DEGs based on p-value and logFC thresholds
  da_results <- da_results |>
    dplyr::filter(!!sym(deg_pval_col) < deg_pval_thresh,
                  abs(!!sym(deg_logfc_col)) > deg_logfc_thresh)

  if (nrow(da_results) == 0) {
    stop("No DEGs meet the specified p-value and logFC thresholds.")
  }

  # Get numeric exposure variables
  numeric_exposures <- expomicset |>
    MultiAssayExperiment::colData() |>
    as.data.frame() |>
    select_if(is.numeric) |>
    colnames()
  if (!is.null(exposure_cols)) {
    numeric_exposures <- intersect(numeric_exposures, exposure_cols)
  }
  if (length(numeric_exposures) == 0) {
    stop("No numeric exposure variables found in colData.")
  }

  # Initialize correlation results storage
  correlation_results <- list()

  # Iterate through each assay
  for (experiment_name in unique(da_results$exp_name)) {
    message("Processing experiment: ", experiment_name)

    if(robust){

      if(!is.null(score_thresh)){

        # Extract relevant DEGs for this assay
        selected_features <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$feature_stability |>
          dplyr::filter(stability_score > score_thresh) |>
          dplyr::filter(exp_name == experiment_name) |>
          dplyr::pull(feature) |>
          unique()

      } else{
        # Use default score_thresh from metadata
        score_thresh <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$score_thresh

        # Extract relevant DEGs for this assay
        selected_features <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$feature_stability |>
          dplyr::filter(stability_score > score_thresh) |>
          dplyr::filter(exp_name == experiment_name) |>
          dplyr::pull(feature) |>
          unique()
      }

    }else{
      # Extract relevant DEGs for this assay
      selected_features <- da_results |>
        dplyr::filter(exp_name == experiment_name) |>
        dplyr::pull(feature) |>
        unique()
    }

    if (length(selected_features) == 0) {
      warning("No relevant DEGs found for ", experiment_name, ", skipping.")
      next
    }

    # Extract and update SummarizedExperiment
    se <- .update_assay_colData(expomicset, experiment_name)

    # Filter SummarizedExperiment to Selected Features
    se <- se[rownames(se) %in% selected_features, , drop = FALSE]

    if (nrow(se) == 0) {
      warning("No valid DEGs left in assay after filtering for ", experiment_name, ", skipping.")
      next
    }

    # Batch Processing
    feature_batches <- split(
      selected_features,
      ceiling(seq_along(selected_features) / batch_size))

    batch_results <- list()

    batch_index <- 1

    for (batch in feature_batches) {
      message(
        "  - Processing batch ",
        batch_index,
        " of ",
        length(feature_batches),
        " (",
        length(batch),
        " features)...")

      batch_index <- batch_index + 1

      # Subset Data to Only Batch Features
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

  # Combine results and adjust p-values
  combined_results <- correlation_results |>
    dplyr::bind_rows() |>
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
    dplyr::left_join(MultiAssayExperiment::metadata(expomicset)$var_info,
                     by=c("exposure"="variable"))


    # dplyr::left_join(expomicset |>
    #                    pivot_feature() |>
    #                    as.data.frame(),
    #                  by=c("feature"=".feature"))

  if (nrow(combined_results) == 0) {
    warning("No significant correlations found in any experiment.")
    return(expomicset)
  }

  if(action=="add"){
    # Save to metadata
    MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation <- combined_results
    return(expomicset)
  }else if (action =="get"){
    return(combined_results)
  }else{
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  message("DEG-exposure correlation analysis completed.")
}
