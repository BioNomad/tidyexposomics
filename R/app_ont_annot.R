#' Build the Ontology Annotation Shiny app
#'
#' Returns a Shiny app object for ontology annotation. This function does not
#' launch the app. Please call \code{shiny::runApp()} yourself (see examples).
#'
#' @param use_demo Logical, if TRUE, load packaged lightweight demo ontologies.
#' @param ... Optional named overrides passed as data.frames with columns
#'   \code{id}, \code{name}, \code{ancestors}, \code{parents}, \code{children}:
#'   \code{hpo}, \code{ecto}, \code{chebi}.
#'
#' @return A \code{shiny.appobj}.
#' @examples
#' if (interactive()) {
#'   app <- build_ont_annot_app()
#'   shiny::runApp(app)
#' }
#' @export
#' @importFrom shiny shinyApp
build_ont_annot_app <- function(use_demo = TRUE, ...) {
  ont <- .load_ontologies(use_demo = use_demo, ...)
  ui  <- .ont_annot_build_ui()
  srv <- .ont_annot_build_server(ontologies = ont)
  shiny::shinyApp(ui = ui, server = srv)
}
