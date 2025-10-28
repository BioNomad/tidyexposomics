#' Summarize and Visualize Analysis Pipeline Steps
#'
#' This function prints and visualizes the analysis steps stored in the
#' metadata of a \code{MultiAssayExperiment} object. The steps are optionally
#' printed to the console as a numbered list and can be rendered as a
#' left-to-right Mermaid flowchart.
#' The flowchart connects steps with arrows and includes step notes
#' if requested.
#'
#' @param exposomicset A \code{MultiAssayExperiment} object that contains a
#' "summary"
#' entry in its metadata, which includes a list named \code{steps}.
#' @param show_index Logical, default \code{TRUE}. If \code{TRUE},
#' prefixes each step with its index.
#' @param console_print Logical, default \code{TRUE}. If \code{TRUE},
#' prints the step list to the console.
#' @param diagram_print Logical, default \code{FALSE}. If \code{TRUE},
#'  renders a Mermaid diagram of the steps.
#' @param include_notes Logical, default \code{TRUE}. If \code{TRUE},
#' appends any "notes" associated with each step to the label.
#'
#' @return No return value. This function is called for its side effects:
#' console output and/or diagram rendering.
#'
#' @details The Mermaid flowchart is rendered left-to-right and connects
#' each step in sequence. Each node is labeled using the step name and,
#'  optionally,
#' any attached notes.
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # Test for normality
#' mae <- mae |>
#'     run_normality_check() |>
#'     transform_exposure(exposure_cols = c("age", "bmi", "exposure_pm25"))
#'
#' # Run the pipeline summary
#' run_pipeline_summary(mae)
#'
#' @export
run_pipeline_summary <- function(
  exposomicset,
  show_index = TRUE,
  console_print = TRUE,
  diagram_print = FALSE,
  include_notes = TRUE
) {
    # Validate metadata
    summary_md <- MultiAssayExperiment::metadata(exposomicset)$summary
    if (!("steps" %in% names(summary_md))) {
        stop("Please run analysis steps before running the pipeline summary.")
    }

    # Get step names
    step_names <- names(summary_md$steps)

    # Optionally include notes in step text
    step_labels <- purrr::map_chr(step_names, function(step) {
        note <- if (include_notes && !is.null(summary_md$steps[[step]]$notes)) {
            paste0(" - ", summary_md$steps[[step]]$notes)
        } else {
            ""
        }
        paste0(step, note)
    })

    # Add index if needed
    labeled_steps <- if (show_index) {
        sprintf("%d. %s", seq_along(step_labels), step_labels)
    } else {
        step_labels
    }

    # Optional console output
    if (console_print) {
        message(paste(labeled_steps, collapse = "\n"))
    }

    if (diagram_print) {
        .check_suggested("DiagrammeR")
        # Build Mermaid diagram
        mermaid_lines <- c("graph TD")
        for (i in seq_along(labeled_steps)[-length(labeled_steps)]) {
            from <- sprintf("step%d", i)
            to <- sprintf("step%d", i + 1)
            from_label <- labeled_steps[i]
            to_label <- labeled_steps[i + 1]

            mermaid_lines <- c(
                mermaid_lines,
                sprintf('%s["%s"] --> %s["%s"]', from, from_label, to, to_label)
            )
        }

        DiagrammeR::mermaid(paste(mermaid_lines, collapse = "\n"))
    }
}
