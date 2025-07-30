#' Filter Non-Normal Exposure Variables
#'
#' Removes exposure variables that deviate significantly from a
#' normal distribution based on
#' normality test results stored in metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object containing
#' exposure and omics data.
#' @param p_thresh A numeric value specifying the p-value
#' threshold for normality. Variables with `p.value < p_thresh` are removed.
#'  Default is `0.05`.
#'
#' @details
#' The function identifies exposure variables that fail a normality test
#'  using `metadata(expomicset)$transformation$norm_df`.
#' - Exposure variables with `p.value < p_thresh` are removed from
#' `colData(expomicset)`.
#' - The same filtering is applied to `colData` in each experiment
#' within `experiments(expomicset)`.
#'
#' @return A `MultiAssayExperiment` object with non-normal
#' exposure variables removed.
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'   n_samples = 20,
#'   return_mae = TRUE
#' )
#'
#' # Test for normality
#' mae <- mae |>
#'   run_normality_check() |>
#'   transform_exposure(exposure_cols=c("age","bmi","exposure_pm25"))
#'
#' # Remove non-normal variables
#' mae_filtered <- mae |>
#'   filter_non_normal()
#'
#' @export
filter_non_normal <- function(expomicset,
                              p_thresh = 0.05) {

  if("transformation" %in% names(expomicset |>
                                 MultiAssayExperiment::metadata() |>
                                 purrr::pluck("quality_control"))){
    norm_df <- expomicset |>
      MultiAssayExperiment::metadata() |>
      purrr::pluck("quality_control",
                   "transformation",
                   "norm_df")
  } else{
    norm_df <- expomicset |>
      MultiAssayExperiment::metadata() |>
      purrr::pluck("quality_control",
                   "normality",
                   "norm_df")
  }
  norm_df <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("quality_control",
                 "transformation",
                 "norm_df")

  if (!all(c("exposure", "p.value") %in% colnames(norm_df))) {
    stop("norm_df must contain 'exposure' and 'p.value' columns.")
  }

  # Get exposures to drop
  non_normal_vars <- norm_df |>
    dplyr::filter(p.value < p_thresh) |>
    dplyr::pull(exposure)

  message("Filtering out ",
          length(non_normal_vars),
          " non-normal exposure variables.")

  # Filter colData in the main object
  MultiAssayExperiment::colData(expomicset) <- MultiAssayExperiment::colData(
    expomicset)[
    , !colnames(MultiAssayExperiment::colData(expomicset)) %in% non_normal_vars,
    drop = FALSE]

  # Filter colData in each experiment
  for (omics_name in names(MultiAssayExperiment::experiments(expomicset))) {
    current_colData <- MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(expomicset)[[omics_name]])

    current_colData <- current_colData[
      , !colnames(current_colData) %in% non_normal_vars,
      drop = FALSE]

    MultiAssayExperiment::colData(
      MultiAssayExperiment::experiments(
        expomicset)[[omics_name]]) <- current_colData
  }

  # Log metadata
  step_record <- list(
    filter_non_normal = list(
      timestamp = Sys.time(),
      params = list(p_thresh = p_thresh),
      notes = paste0(
        "Filtered out non-normal exposure variables,
        Based on p-value threshold of ",
        p_thresh, ". ",
        length(non_normal_vars),
        " variables removed.")
    ))

  MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$summary$steps,
    step_record
  )

  return(expomicset)
}

