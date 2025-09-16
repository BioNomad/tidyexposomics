#' Summarize Exposure Variables
#'
#' Computes summary statistics for numeric exposure variables and
#' optionally stores the results in the `MultiAssayExperiment` metadata.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing exposure
#' data in the sample metadata.
#' @param exposure_cols A character vector of exposure variable names
#' to summarize. If `NULL`, all numeric exposure variables are included.
#' @param action A string specifying the action to take. Use `"add"`
#' to attach the summary table to `metadata(exposomicset)` or `"get"`
#' to return the summary table directly. Default is `"add"`.
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
#' - Merges the result with variable metadata stored in
#' `metadata(exposomicset)$codebook`.
#'
#' @return
#' A modified `MultiAssayExperiment` object (if `action = "add"`),
#' or a data frame of summary statistics (if `action = "get"`).
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # Summarize exposure data
#' exp_sum <- mae |>
#'     run_summarize_exposures(
#'         exposure_cols = c("age", "bmi", "exposure_pm25"),
#'         action = "get"
#'     )
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom dplyr select all_of where group_by summarise ungroup
#' inner_join mutate_if everything
#' @importFrom tidyr pivot_longer
#' @importFrom purrr pluck
#' @export
run_summarize_exposures <- function(
    exposomicset,
    exposure_cols = NULL,
    action = "add") {
    # library(dplyr)
    # library(tidyr)

    # Extract exposure data
    exposure_df <- exposomicset |>
        pivot_sample()

    # Subset to selected columns if specified
    if (!is.null(exposure_cols)) {
        exposure_df <- exposure_df |>
            dplyr::select(dplyr::all_of(exposure_cols))
    }

    # Keep only numeric columns
    exposure_df <- exposure_df |>
        dplyr::select(dplyr::where(is.numeric))

    # Pivot to long format for summary
    long_df <- exposure_df |>
        tidyr::pivot_longer(
            cols = dplyr::everything(),
            names_to = "variable",
            values_to = "value"
        )

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
            exposomicset |>
                MultiAssayExperiment::metadata() |>
                purrr::pluck("codebook"),
            by = "variable"
        ) |>
        dplyr::mutate_if(is.numeric, ~ round(., digits = 2))

    if (action == "add") {
        # Store results
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$quality_control$exposure_summary_df <- exposure_summary_df
        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

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
                        paste0(
                            "Subset to user-specified exposure columns (n = ",
                            length(exposure_cols),
                            ")."
                        )
                    }
                )
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else if (action == "get") {
        # Return the summary data frame
        return(exposure_summary_df)
    } else {
        stop("Invalid action. Use 'add' or 'get'.")
    }
}
