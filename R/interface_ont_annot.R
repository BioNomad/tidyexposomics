#' Internal - builds ontology annotation UI
#'
#' @keywords internal
#' @return A \code{shiny.tag} UI definition.
#' @importFrom shiny fluidPage tags HTML div img titlePanel fileInput uiOutput
.ont_annot_build_ui <- function() {
    logo_path <- system.file("app", "www", "logo.png", package = "tidyexposomics")
    shiny::fluidPage(
        shiny::tags$head(shiny::tags$style(shiny::HTML(
            "body { background:#fff; color:#1e0c1f; font-family:Helvetica,Arial,sans-serif; }
       h1,h4 { color:#1e0c1f; }
       .btn { background:#71255A; border-color:#71255A; color:#fff; }
       .btn:hover { background:#5a1e47; border-color:#5a1e47; }
       .selectize-input, .form-control { background:#fff; color:#1e0c1f; }
       .dataTables_wrapper .paginate_button { color:#71255A!important; }
       .shiny-input-container { margin-bottom:20px; }"
        ))),
        shiny::div(
            style = "position:absolute;top:10px;right:10px;z-index:1000;",
            if (file.exists(logo_path)) shiny::img(src = "www/logo.png", height = "150px")
        ),
        shiny::titlePanel(shiny::div("Ontology Annotation Tool", style = "color:#1e0c1f")),
        shiny::fileInput("user_df", "Upload Your Exposure Metadata (CSV)", accept = ".csv"),
        shiny::uiOutput("main_ui")
    )
}
