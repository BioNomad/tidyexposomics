#' Filter Features and Variables with High Missingness
#'
#' Removes exposure variables and omics features with missing values above a specified threshold.
#' Generates missing data summaries and quality control (QC) plots.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure and omics data.
#' @param na_thresh A numeric value specifying the percentage of missing data allowed before a variable or feature is removed. Default is `20`.
#' @param na_plot_thresh A numeric value specifying the minimum missing percentage for inclusion in QC plots. Default is `5`.
#'
#' @details
#' The function assesses missingness in both `colData(expomicset)` (exposure data) and `experiments(expomicset)` (omics data).
#' - Exposure variables with more than `na_thresh`% missing values are removed.
#' - Omics features (rows in assay matrices) exceeding `na_thresh`% missing values are filtered.
#' - Missingness summaries and QC plots are generated using `naniar::gg_miss_var()` and stored in metadata.
#'
#' @return A `MultiAssayExperiment` object with filtered exposure variables and omics features.
#' QC results, including missingness summaries and plots, are stored in `metadata(expomicset)$na_qc`.
#'
#' @examples
#' \dontrun{
#' filtered_expom <- filter_missing(
#'   expomicset = expom,
#'   na_thresh = 20,
#'   na_plot_thresh = 5
#' )
#' }
#'
#' @export
filter_missing <- function(expomicset, 
                           na_thresh = 20, 
                           na_plot_thresh = 5) {
  
  # Helper function to identify variables/features to exclude
  get_vars_to_exclude <- function(data, thresh) {
    miss_summary <- naniar::miss_var_summary(data)
    list(
      to_exclude = miss_summary |>
        dplyr::filter(pct_miss > thresh) |> 
        dplyr::pull(variable),
      summary = miss_summary
    )
  }
  
  # Handle metadata (colData)
  exposure <- as.data.frame(MultiAssayExperiment::colData(expomicset))
  col_exclude <- get_vars_to_exclude(exposure, na_thresh)
  MultiAssayExperiment::colData(expomicset) <- as(
    MultiAssayExperiment::colData(expomicset)[
      , !colnames(MultiAssayExperiment::colData(expomicset)) %in% 
        col_exclude$to_exclude,
      drop = FALSE], 
    "DataFrame")
  
  # QC plot for metadata
  MultiAssayExperiment::metadata(expomicset)$na_qc <- list(
    exposure = list(
      vars_to_exclude_exposure_sum = col_exclude$summary |> 
        dplyr::filter(pct_miss > na_thresh),
      all_var_sum = col_exclude$summary |> 
        dplyr::filter(pct_miss > 0),
      na_exposure_qc_plot = naniar::gg_miss_var(
        exposure |> 
          dplyr::select(
            dplyr::all_of(col_exclude$to_exclude))
      )
    )
  )
  
  message("Missing Data Filter threshold: ", na_thresh, "%")
  message("Filtered metadata variables: ", paste(col_exclude$to_exclude, collapse = ", "))
  
  # Handle omics experiments
  for (omics_name in names(MultiAssayExperiment::experiments(expomicset))) {
    experiment <- MultiAssayExperiment::experiments(expomicset)[[omics_name]]
    original_assay <- SummarizedExperiment::assays(experiment)[[1]]
    
    # Transpose and convert assay to data.frame for processing
    assay_data <- as.data.frame(t(original_assay))
    
    # Identify features (rows) to exclude based on missingness
    row_exclude <- get_vars_to_exclude(assay_data, na_thresh)
    
    # Filter rows with high missingness
    experiment <- experiment[!rownames(experiment) %in% row_exclude$to_exclude,]
    
    MultiAssayExperiment::experiments(expomicset)[[omics_name]] <- experiment
    
    # QC plot for omics data
    MultiAssayExperiment::metadata(expomicset)$na_qc[[omics_name]] <- list(
      vars_to_exclude_omics_sum = row_exclude$summary |> 
        dplyr::filter(pct_miss > na_thresh),
      all_var_sum = row_exclude$summary |> 
        dplyr::filter(pct_miss > 0),
      na_omics_qc_plot = naniar::gg_miss_var(
        assay_data |> 
          dplyr::select(all_of(row_exclude$to_exclude)))
    )
    
    message("Filtered rows with high missingness in ", omics_name, ": ", length(row_exclude$to_exclude))
  }
  
  return(expomicset)
}

