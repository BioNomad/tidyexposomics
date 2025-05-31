#' Summarize and Visualize Analysis Pipeline Steps
#'
#' This function prints and visualizes the analysis steps stored in the
#' metadata of a `MultiAssayExperiment` object. The steps are optionally printed
#' to the console as a numbered list and always rendered as a left-to-right Mermaid
#' flowchart. The flowchart wraps to a new row after every 5 steps.
#'
#' @param expomicset A `MultiAssayExperiment` object that contains a `"steps"`
#' entry in its metadata. This entry should be a character vector describing the pipeline steps.
#' @param show_index Logical, default `TRUE`. If `TRUE`, each step is prefixed with its index
#' (e.g., "1. Step description"). If `FALSE`, only the step text is shown.
#' @param console_print Logical, default `FALSE`. If `TRUE`, prints the step list to the console.
#'
#' @return No return value. This function is called for its side effects:
#' it optionally prints steps to the console and renders a Mermaid diagram.
#'
#' @examples
#' if (requireNamespace("MultiAssayExperiment") && requireNamespace("DiagrammeR")) {
#'   mae <- MultiAssayExperiment::MultiAssayExperiment()
#'   MultiAssayExperiment::metadata(mae)$steps <- c(
#'     "Import data", "Normalize", "Run PCA", "Fit model", "Export results",
#'     "Visualize", "Report"
#'   )
#'   run_pipeline_summary(mae, show_index = TRUE, console_print = TRUE)
#' }
#'
#' @export
run_pipeline_summary <- function(expomicset, show_index = TRUE, console_print = FALSE) {
  # Check for 'steps' metadata
  if (!("steps" %in% names(MultiAssayExperiment::metadata(expomicset)))) {
    stop("Please run analysis steps before running the pipeline summary.")
  }

  steps <- MultiAssayExperiment::metadata(expomicset)[["steps"]]


  # Optionally add line numbers
  labeled_steps <- if (show_index) {
    sprintf("%d. %s", seq_along(steps), steps)
  } else {
    steps
  }

  if(console_print) {
    # Print steps to console
    cat(labeled_steps, sep = "\n")
  }

  # Build Mermaid flowchart
  mermaid_lines <- c("graph TD")
  for (i in seq_along(steps)[-length(steps)]) {
    from <- sprintf("step%d", i)
    to <- sprintf("step%d", i + 1)
    from_label <- labeled_steps[i]
    to_label <- labeled_steps[i + 1]

    mermaid_lines <- c(
      mermaid_lines,
      sprintf('%s["%s"] --> %s["%s"]', from, from_label, to, to_label)
    )
  }

  # Render Mermaid diagram
  DiagrammeR::mermaid(paste(mermaid_lines, collapse = "\n"))
}

