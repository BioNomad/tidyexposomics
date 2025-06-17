#' Summarize Exposure Variables
#'
#' Computes summary statistics for numeric exposure variables and optionally stores the results in the `MultiAssayExperiment` metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure data in the sample metadata.
#' @param exposure_cols A character vector of exposure variable names to summarize. If `NULL`, all numeric exposure variables are included.
#' @param action A string specifying the action to take. Use `"add"` to attach the summary table to `metadata(expomicset)` or `"get"` to return the summary table directly. Default is `"add"`.
#'
#' @details
#' This function:
#' - Extracts sample-level exposure data using `pivot_sample()`.
#' - Filters to user-specified exposures (`exposure_cols`) if provided.
#' - Computes descriptive statistics for each numeric variable:
#'   - Number of values (`n_values`)
#'   - Number of NAs (`n_na`)
#'   - Minimum, maximum, and range
#'   - Sum, median, mean
#'   - Standard error of the mean
#'   - 95% confidence interval of the mean
#'   - Variance, standard deviation
#'   - Coefficient of variation (`sd / mean`)
#' - Merges the result with variable metadata stored in `metadata(expomicset)$var_info`.
#'
#' @return
#' A modified `MultiAssayExperiment` object (if `action = "add"`), or a data frame of summary statistics (if `action = "get"`).
#'
#' @examples
#' \dontrun{
#' # Add summaries to metadata
#' expom <- run_summarize_exposures(expom)
#'
#' # Retrieve just the summary table
#' summary_df <- run_summarize_exposures(expom, action = "get")
#' }
#'
#' @export

run_summarize_exposures <- function(
    expomicset,
    exposure_cols = NULL,
    action = "add"
){
  library(dplyr)
  library(tidyr)

  # Extract exposure data
  exposure_df <- expomicset |>
    pivot_sample()

  # Subset to selected columns if specified
  if(!is.null(exposure_cols)){
    exposure_df <- exposure_df |>
      dplyr::select(dplyr::all_of(exposure_cols))
  }

  # Keep only numeric columns
  exposure_df <- exposure_df |>
    dplyr::select(dplyr::where(is.numeric))

  # Pivot to long format for summary
  long_df <- exposure_df |>
    tidyr::pivot_longer(cols = dplyr::everything(),
                 names_to = "variable",
                 values_to = "value")

  # Group by variable and compute stats
  exposure_summary_df <- long_df |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n_values = sum(!is.na(value)),
      n_na = sum(is.na(value)),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      range = max - min,
      sum = sum(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(n_values),
      ci_lower = mean - qt(0.975, df = n_values - 1) * se,
      ci_upper = mean + qt(0.975, df = n_values - 1) * se,
      variance = var(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      coef_var = sd / mean
    ) |>
    dplyr::ungroup()

  # Merge with variable info
  exposure_summary_df <- exposure_summary_df |>
    dplyr::inner_join(
      expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("var_info"),
    by="variable") |>
    dplyr::mutate_if(is.numeric,~round(.,digits = 2))

  if(action=="add"){
    # Store results
    MultiAssayExperiment::metadata(expomicset)$quality_control$exposure_summary_df <- exposure_summary_df

    # Add step record
    step_record <- list(
      run_summarize_exposures = list(
        timestamp = Sys.time(),
        params = list(
          exposure_cols = exposure_cols,
          n_vars_summarized = length(unique(exposure_summary_df$variable))
        ),
        notes = paste0(
          "Summarized ", length(unique(exposure_summary_df$variable)),
          " numeric exposure variables. ",
          if (is.null(exposure_cols)) {
            "Included all numeric exposures from colData."
          } else {
            paste0("Subset to user-specified exposure columns (n = ", length(exposure_cols), ").")
          }
        )
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)

  }else if (action=="get"){
    # Return the summary data frame
    return(exposure_summary_df)

  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
