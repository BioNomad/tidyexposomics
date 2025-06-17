#' Summarize and Visualize Analysis Pipeline Steps
#'
#' This function prints and visualizes the analysis steps stored in the
#' metadata of a `MultiAssayExperiment` object. The steps are optionally printed
#' to the console as a numbered list and always rendered as a left-to-right Mermaid
#' flowchart. The flowchart wraps to a new row after every 5 steps.
#'
#' @param expomicset A `MultiAssayExperiment` object that contains a `"steps"`
#' entry in its metadata.
#' @param show_index Logical, default `TRUE`. Prefix steps with index.
#' @param console_print Logical, default `FALSE`. If `TRUE`, prints to console.
#' @param include_notes Logical, default `FALSE`. If `TRUE`, appends notes from each step.
#'
#' @return No return value; called for side effects (console print + diagram).
#' @export
run_pipeline_summary <- function(expomicset,
                                 show_index = TRUE,
                                 console_print = FALSE,
                                 include_notes = FALSE) {
  # Validate metadata
  summary_md <- MultiAssayExperiment::metadata(expomicset)$summary
  if (!("steps" %in% names(summary_md))) {
    stop("Please run analysis steps before running the pipeline summary.")
  }

  # Get step names
  step_names <- names(summary_md$steps)

  # Optionally include notes in step text
  step_labels <- purrr::map_chr(step_names, function(step) {
    note <- if (include_notes && !is.null(summary_md$steps[[step]]$notes)) {
      paste0(" — ", summary_md$steps[[step]]$notes)
    } else ""
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
    cat(labeled_steps, sep = "\n")
  }

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
