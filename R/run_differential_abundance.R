#' Run Differential Abundance Analysis
#'
#' Performs differential abundance testing across all assays in a
#' `MultiAssayExperiment` object using a specified statistical method
#' such as `limma_voom`. The function updates each assay with its
#' corresponding `colData`, fits the model using the provided formula,
#' and combines the results into a unified table.
#'
#' @param expomicset A `MultiAssayExperiment` containing assays to test.
#' @param formula A model formula for the differential analysis
#' (e.g., ~ group + batch).
#' @param abundance_col Character. The name of the assay matrix to use
#' for abundance values. Default is `"counts"`.
#' @param method Character. Differential analysis method to use.
#' Currently supports `"limma_voom"` (default).
#' @param contrasts A named list of contrasts for pairwise comparisons.
#'  Default is `NULL` (uses default group comparisons).
#' @param scaling_method Character. Scaling method to apply before modeling.
#' Options include `"none"` (default), `"zscore"`, etc.
#' @param action Character. Whether to `"add"` results to `expomicset`
#'  metadata or `"get"` the results as a data frame. Default is `"add"`.
#'
#' @return Either the updated `MultiAssayExperiment` (if `action = "add"`)
#'  or a tibble with differential abundance results (if `action = "get"`).
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 10,
#'   return_mae=TRUE
#'   )
#'
#' # perform differential abundance analysis
#' mae <- run_differential_abundance(
#'   expomicset = mae,
#'   formula = ~ smoker + sex,
#'   abundance_col = "counts",
#'   method = "limma_voom",
#'   action = "add"
#' )
#'
#' @export
run_differential_abundance <- function(
    expomicset,
    formula,
    abundance_col = "counts",
    method = "limma_voom",
    contrasts = NULL,
    scaling_method = "none",
    action="add"
) {

  message("Running differential abundance testing.")

  # Initialize a data frame to store results
  da_results_df <- list()

  # Iterate through assays in expomicset
  for (exp_name in names(MultiAssayExperiment::experiments(expomicset))) {
    message("Processing assay: ", exp_name)

    # Update assay with colData
    exp <- .update_assay_colData(expomicset, exp_name)

    # Run differential analysis using `.run_se_differential_abundance`
    res <- .run_se_differential_abundance(
      se = exp,
      formula = formula,
      abundance_col = abundance_col,
      method = method,
      scaling_method = scaling_method,
      contrasts = contrasts
    )

    # If results exist, append assay name and store them
    if (!is.null(res) && nrow(res) > 0) {
      res <- res |>
        dplyr::mutate(exp_name = exp_name)
      da_results_df[[exp_name]] <- res
    } else {
      warning("No significant results found for assay: ", exp_name)
    }
  }

  # Combine results across assays
  final_results <- da_results_df |>
    dplyr::bind_rows() |>
    as_tibble()

  message("Differential abundance testing completed.")

  if(action == "add") {
    all_metadata <- MultiAssayExperiment::metadata(expomicset)
    all_metadata$differential_analysis$differential_abundance <- final_results
    MultiAssayExperiment::metadata(expomicset) <- all_metadata

    # Add step record to metadata
    step_record <- list(
      run_differential_abundance = list(
        timestamp = Sys.time(),
        params = list(
          formula = format(formula),
          method = method,
          scaling_method = scaling_method,
          abundance_col = abundance_col,
          contrasts = contrasts
        ),
        notes = "Performed differential abundance analysis across all assays."
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )
    return(expomicset)


  } else if (action == "get") {
    return(final_results)
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
  return(expomicset)
}
