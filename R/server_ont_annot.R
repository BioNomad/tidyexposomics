#' Internal - builds ontology annotation server
#'
#' @keywords internal
#' @param ontologies A list with \code{hpo}, \code{ecto}, \code{chebi} data.frames.
#' @return A Shiny server function.
#' @importFrom shiny renderUI column fluidRow h4 selectizeInput actionButton hr
#' @importFrom shiny selectInput numericInput downloadButton observeEvent
#' @importFrom shiny reactiveValues showModal modalDialog updateSelectizeInput
#' @importFrom shiny req HTML
#' @importFrom DT datatable renderDT DTOutput dataTableProxy replaceData formatStyle styleEqual
#' @importFrom dplyr mutate distinct bind_rows case_when filter any_of
.ont_annot_build_server <- function(ontologies) {
    force(ontologies)

    preprocess_ont <- function(ont) {
        fix_rel <- function(x) if (!is.list(x)) strsplit(x, ";\\s*") else x
        dplyr::mutate(
            ont,
            ancestors = fix_rel(ancestors),
            parents   = fix_rel(parents),
            children  = fix_rel(children),
            depth     = lengths(ancestors)
        )
    }

    # Prepare once
    ont_list <- lapply(ontologies[c("hpo", "ecto", "chebi")], preprocess_ont)

    ontology_df_all <-
        dplyr::bind_rows(
            dplyr::mutate(ont_list$hpo[, c("id", "name")], source = "hpo"),
            dplyr::mutate(ont_list$ecto[, c("id", "name")], source = "ecto"),
            dplyr::mutate(ont_list$chebi[, c("id", "name")], source = "chebi")
        ) |>
        dplyr::distinct()

    function(input, output, session) {
        rv <- shiny::reactiveValues(data = NULL)

        shiny::observeEvent(input$user_df, {
            shiny::req(input$user_df)
            df <- utils::read.csv(input$user_df$datapath, stringsAsFactors = FALSE)

            if (!"variable" %in% names(df)) {
                shiny::showModal(shiny::modalDialog("Your file needs a 'variable' column.", easyClose = TRUE))
                return()
            }

            expected <- c(
                "selected_ontology_label", "selected_ontology_id",
                "root_id", "root_label", "category", "category_source"
            )
            for (col in expected) if (!col %in% names(df)) df[[col]] <- NA_character_

            rv$data <- df
        })

        output$main_ui <- shiny::renderUI({
            shiny::req(rv$data)
            shiny::fluidRow(
                shiny::column(
                    4,
                    shiny::h4("Step 1: Annotate Variable"),
                    shiny::uiOutput("selected_var"),
                    shiny::selectizeInput("ontology_choice", "Choose Term:",
                        choices = NULL,
                        options = list(placeholder = "Type to search")
                    ),
                    shiny::actionButton("apply_annotation", "Apply Annotation"),
                    shiny::uiOutput("ontology_description"),
                    shiny::hr(),
                    shiny::h4("Step 2: Categorize Variables"),
                    shiny::selectInput("cat_ontology", "Ontology:", choices = c("hpo", "ecto", "chebi"), selected = "hpo"),
                    shiny::numericInput("cat_root_depth", "Root depth level:", value = 2, min = 0),
                    shiny::actionButton("apply_categorization", "Apply Categorization"),
                    shiny::hr(),
                    shiny::h4("Step 3: Download Annotated Data"),
                    shiny::hr(),
                    shiny::h4("Legend"),
                    shiny::div(
                        style = "padding-left:10px;",
                        shiny::tags$div(
                            style = "background:#fde0dc;padding:4px;border-radius:3px;margin-bottom:4px;",
                            "Missing ontology annotation"
                        ),
                        shiny::tags$div(
                            style = "background:#fff2cc;padding:4px;border-radius:3px;margin-bottom:4px;",
                            "Missing category"
                        ),
                        shiny::tags$div(
                            style = "background:#fff7e6;padding:4px;border-radius:3px;",
                            "Manually categorized"
                        )
                    ),
                    shiny::downloadButton("save_csv", "Download Annotated CSV")
                ),
                shiny::column(8, DT::DTOutput("table"))
            )
        })

        shiny::observe({
            shiny::req(rv$data)
            shiny::updateSelectizeInput(session, "ontology_choice",
                choices = ontology_df_all$name, server = TRUE
            )
        })

        output$selected_var <- shiny::renderUI({
            sel <- input$table_rows_selected
            if (length(sel)) {
                shiny::HTML(paste0("<strong>Selected:</strong> ", rv$data$variable[sel]))
            } else {
                shiny::HTML("<em>No row selected</em>")
            }
        })

        output$table <- DT::renderDT(
            {
                shiny::req(rv$data)
                display_data <- rv$data
                display_data$row_status <- dplyr::case_when(
                    is.na(display_data$selected_ontology_id) ~ "Missing Annotation",
                    is.na(display_data$category) ~ "Missing Category",
                    display_data$category_source == "manual" ~ "Manual",
                    TRUE ~ "OK"
                )

                DT::datatable(
                    display_data,
                    selection = "multiple",
                    rownames = FALSE,
                    editable = list(
                        target = "cell",
                        disable = list(columns = which(names(display_data) != "category") - 1)
                    ),
                    options = list(
                        pageLength = 2000,
                        scrollY    = "400px",
                        scrollX    = TRUE,
                        autoWidth  = TRUE,
                        stateSave  = TRUE,
                        columnDefs = list(list(visible = FALSE, targets = which(names(display_data) == "row_status") - 1))
                    ),
                    class = "stripe nowrap"
                ) |>
                    DT::formatStyle("row_status",
                        target = "row",
                        backgroundColor = DT::styleEqual(
                            c("Missing Annotation", "Missing Category", "Manual"),
                            c("#fde0dc", "#fff2cc", "#fff7e6")
                        )
                    )
            },
            server = TRUE
        )

        proxy <- DT::dataTableProxy("table")

        shiny::observeEvent(input$apply_annotation, {
            shiny::req(rv$data, input$ontology_choice)
            rows <- input$table_rows_selected
            if (length(rows)) {
                sel <- dplyr::filter(ontology_df_all, name == input$ontology_choice)
                if (nrow(sel) == 1) {
                    rv$data$selected_ontology_label[rows] <- sel$name
                    rv$data$selected_ontology_id[rows] <- sel$id
                }
            }
            DT::replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
        })

        output$ontology_description <- shiny::renderUI({
            shiny::req(input$ontology_choice)
            term <- dplyr::filter(ontology_df_all, name == input$ontology_choice)
            if (nrow(term) == 1) {
                desc <- .get_ols_description(term$id, sub(":.*$", "", term$id))
                if (!is.na(desc)) {
                    shiny::div(
                        style = "border-left:4px solid #71255A;padding-left:10px;margin:10px 0;",
                        shiny::strong("Description:"), shiny::br(), desc
                    )
                }
            }
        })

        shiny::observeEvent(input$apply_categorization, {
            shiny::req(rv$data)
            if (all(is.na(rv$data$selected_ontology_id))) {
                shiny::showModal(shiny::modalDialog("Annotate at least one ID before categorizing.", easyClose = TRUE))
                return()
            }
            rows <- input$table_rows_selected
            if (!length(rows)) {
                shiny::showModal(shiny::modalDialog("Please select one or more rows to categorize.", easyClose = TRUE))
                return()
            }

            onto_name <- input$cat_ontology
            root_lvl <- as.numeric(input$cat_root_depth)
            target_df <- ont_list[[onto_name]]

            updated <- .run_categorize_ontology(
                data         = rv$data[rows, , drop = FALSE],
                id_col       = "selected_ontology_id",
                ontologyDF   = target_df,
                root_level   = root_lvl,
                assign_label = TRUE
            )

            rv$data$root_id[rows] <- updated$root_id
            rv$data$root_label[rows] <- updated$root_label
            rv$data$category[rows] <- updated$category
            rv$data$category_source[rows] <- updated$category_source

            DT::replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
        })

        shiny::observeEvent(input$table_cell_edit, {
            info <- input$table_cell_edit
            i <- info$row
            j <- info$col
            v <- info$value
            colname <- names(rv$data)[j + 1] # DT is 0-based

            if (identical(colname, "category") &&
                (is.na(rv$data$root_id[i]) || nchar(rv$data$root_id[i]) <= 1)) {
                rv$data$category[i] <- v
                rv$data$category_source[i] <- "manual"
                DT::replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
            } else {
                shiny::showModal(shiny::modalDialog(
                    title = "Editing not allowed",
                    "Only 'Unmapped' categories can be changed manually.",
                    easyClose = TRUE
                ))
            }
        })

        output$save_csv <- shiny::downloadHandler(
            filename = function() paste0("annotated_", Sys.Date(), ".csv"),
            content  = function(file) utils::write.csv(rv$data, file, row.names = FALSE)
        )
    }
}
