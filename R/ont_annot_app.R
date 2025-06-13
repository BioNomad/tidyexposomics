#' Launch Ontology Annotation Shiny App
#'
#' This function launches the Shiny application for annotating exposome variable names
#' using ontology terms such as HPO, ENVO, and ChEBI. The app is bundled with the
#' `tidyexposomics` package and allows interactive selection or search of ontology terms.
#'
#' @return Launches a Shiny app in a new window; does not return a value.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   ont_annot_app()
#' }
ont_annot_app <- function() {
  app_dir <- system.file("shiny/ont_annot", package = "tidyexposomics")
  shiny::runApp(app_dir, display.mode = "normal")
}
