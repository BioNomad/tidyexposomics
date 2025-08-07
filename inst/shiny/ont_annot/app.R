library(shiny)
library(tidyexposomics)
library(DT)
library(dplyr)
library(ontologyIndex)
library(httr)
library(jsonlite)

# ---- Load lightweight annotation ontology ----
# data("hpo", package = "tidyexposomics")
# data("ecto", package = "tidyexposomics")
# data("chebi", package = "tidyexposomics")

load_annotation_data()

#data("ont")
# ontology_df <- ontology_df |>
#   dplyr::select(id, name) |>
#   dplyr::distinct()

# ---- Load & preprocess the full ontologies ONCE ----
preprocess_ont <- function(ont) {
  fix_rel <- function(x) if (!is.list(x)) strsplit(x, ";\\s*") else x
  ont %>%
    mutate(
      ancestors = fix_rel(ancestors),
      parents   = fix_rel(parents),
      children  = fix_rel(children),
      depth     = lengths(ancestors)
    )
}

# hpo_raw   <- readRDS("../../../data/hpo.rds")
# ecto_raw  <- readRDS("../../../data/ecto.rds")
# chebi_raw <- readRDS("../../../data/chebi.rds")

# hpo_raw   <- data("hpo")
# ecto_raw  <- data("ecto")
# chebi_raw <- data("chebi")

# hpo_raw   <- get("hpo", envir = asNamespace("tidyexposomics"))
# ecto_raw  <- get("ecto", envir = asNamespace("tidyexposomics"))
# chebi_raw <- get("chebi", envir = asNamespace("tidyexposomics"))

hpo_raw   <- hpo
ecto_raw  <- ecto
chebi_raw <- chebi

ont_list <- list(
  hpo   = preprocess_ont(hpo_raw),
  ecto  = preprocess_ont(ecto_raw),
  chebi = preprocess_ont(chebi_raw)
)

# ---- Helper: Fetch OLS description ----
get_ols_description <- function(ontology_id, ontology_prefix) {
  iri <- paste0("http://purl.obolibrary.org/obo/", gsub(":", "_", ontology_id))
  double_encode <- function(url) {
    tmp <- URLencode(url, reserved = TRUE)
    gsub("%", "%25", tmp, fixed = TRUE)
  }
  enc <- double_encode(iri)
  url <- paste0("https://www.ebi.ac.uk/ols4/api/ontologies/",
                tolower(ontology_prefix), "/terms/", enc)
  res <- try(GET(url), silent = TRUE)
  if (inherits(res, "try-error") || status_code(res) != 200) return(NA_character_)
  j <- fromJSON(rawToChar(res$content))
  desc <- j$description
  if (!is.null(desc) && length(desc)>0 && nzchar(desc[[1]])) {
    desc[[1]]
  } else if (!is.null(j$annotation$description) && length(j$annotation$description)>0) {
    j$annotation$description[[1]]
  } else {
    "No description found from OLS."
  }
}

# ---- Categorization function (uses ont_list) ----
run_categorize_ontology <- function(
    data,
    id_col,
    ontology,
    root_level = 0,
    assign_label = TRUE
) {
  ontology_df <- switch(
    ontology,
    "hpo"   = hpo,
    "ecto"  = ecto,
    "chebi" = chebi,
    stop("Invalid ontology. Choose from 'hpo', 'ecto', or 'chebi'.")
  )
  fix_rel_cols <- function(x) if (!is.list(x)) strsplit(x, ";\\s*") else x
  ontology_df <- ontology_df %>%
    dplyr::mutate(
      ancestors = fix_rel_cols(ancestors),
      parents   = fix_rel_cols(parents),
      children  = fix_rel_cols(children),
      depth     = lengths(ancestors)
    )
  id_to_ancestors <- setNames(ontology_df$ancestors, ontology_df$id)
  id_to_depth     <- setNames(ontology_df$depth,     ontology_df$id)
  id_to_name      <- setNames(ontology_df$name,      ontology_df$id)

  assign_to_root <- switch(
    class(root_level),

    # character: explicit IDs
    character = {
      root_nodes <- root_level
      function(term_id) {
        anc <- c(term_id, id_to_ancestors[[term_id]])
        matched <- intersect(anc, root_nodes)
        if (length(matched)) matched[1] else NA_character_
      }
    },

    # numeric: exact depth match
    numeric = {
      function(term_id) {
        anc    <- c(term_id, id_to_ancestors[[term_id]])
        depths <- id_to_depth[anc]
        idxs   <- which(depths == root_level)
        if (length(idxs)) {
          anc[idxs[1]]
        } else {
          NA_character_
        }
      }
    },

    # default: ontology roots (no parents)
    {
      root_nodes <- ontology_df$id[lengths(ontology_df$parents) == 0]
      function(term_id) {
        anc <- c(term_id, id_to_ancestors[[term_id]])
        matched <- intersect(anc, root_nodes)
        if (length(matched)) matched[1] else NA_character_
      }
    }
  )

  term_ids     <- data[[id_col]]
  assigned_ids <- vapply(term_ids, assign_to_root, character(1))
  assigned_labels <- id_to_name[assigned_ids]

  if (assign_label) {
    data %>%
      dplyr::select(-any_of(c("root_id", "root_label", "category", "category_source"))) %>%
      dplyr::mutate(
        root_id         = assigned_ids,
        root_label      = assigned_labels,
        category        = ifelse(is.na(assigned_labels), "Unmapped", assigned_labels),
        category_source = ifelse(is.na(assigned_labels), "manual", "ontology")
      )
  } else {
    assigned_ids
  }
}

# ---- UI ----
ui <- fluidPage(
  tags$head(tags$style(HTML(
    "body { background:#fff; color:#1e0c1f; font-family:Helvetica,Arial,sans-serif; }
    h1,h4 { color:#1e0c1f; }
    .btn { background:#71255A; border-color:#71255A; color:#fff; }
    .btn:hover { background:#5a1e47; border-color:#5a1e47; }
    .selectize-input, .form-control { background:#fff; color:#1e0c1f; }
    .dataTables_wrapper .paginate_button { color:#71255A!important; }
    .shiny-input-container { margin-bottom:20px; }"
  ))),
  div(style="position:absolute;top:10px;right:10px;z-index:1000;",
      img(src="logo.png", height="150px")),
  titlePanel(div("Ontology Annotation Tool", style="color:#1e0c1f")),
  fileInput("user_df", "Upload Your Exposure Metadata (CSV)", accept=".csv"),
  uiOutput("main_ui")
)

# ---- Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(data = NULL)

  observeEvent(input$user_df, {
    req(input$user_df)
    # df <- read.csv(input$user_df$datapath, stringsAsFactors=FALSE)
    # if (!"variable" %in% names(df)) {
    #   showModal(modalDialog("Your file needs a 'variable' column.", easyClose=TRUE))
    #   return()
    # }
    # df$selected_ontology_label <- df$selected_ontology_label %||% NA_character_
    # df$selected_ontology_id    <- df$selected_ontology_id    %||% NA_character_
    # rv$data <- df

    df <- read.csv(input$user_df$datapath, stringsAsFactors = FALSE)

    # Must include a 'variable' column
    if (!"variable" %in% names(df)) {
      showModal(modalDialog("Your file needs a 'variable' column.", easyClose = TRUE))
      return()
    }

    # Add any missing expected columns
    expected_cols <- c("selected_ontology_label", "selected_ontology_id", "root_id", "root_label", "category", "category_source")
    for (col in expected_cols) {
      if (!col %in% names(df)) {
        df[[col]] <- NA_character_
      }
    }

    rv$data <- df

  })

  output$main_ui <- renderUI({
    req(rv$data)
    fluidRow(
      column(4,
             h4("Step 1: Annotate Variable"),
             uiOutput("selected_var"),
             selectizeInput("ontology_choice","Choose Term:", choices=NULL, options=list(placeholder="Type to search")),
             actionButton("apply_annotation","Apply Annotation"),
             uiOutput("ontology_description"),
             hr(),
             h4("Step 2: Categorize Variables"),
             selectInput("cat_ontology","Ontology:",choices=c("hpo","ecto","chebi"),selected="hpo"),
             numericInput("cat_root_depth","Root depth level:",value=2,min=0),
             actionButton("apply_categorization","Apply Categorization"),
             hr(),
             h4("Step 3: Download Annotated Data"),
             hr(),
             h4("Legend"),
             div(style="padding-left:10px;",
                 tags$div(style="background:#fde0dc;padding:4px;border-radius:3px;margin-bottom:4px;",
                          "Missing ontology annotation"),
                 tags$div(style="background:#fff2cc;padding:4px;border-radius:3px;margin-bottom:4px;",
                          "Missing category"),
                 tags$div(style="background:#fff7e6;padding:4px;border-radius:3px;",
                          "Manually categorized")
             ),
             downloadButton("save_csv","Download Annotated CSV")
      ),
      column(8,
             DTOutput("table")
      )
    )
  })

  observe({
    req(rv$data)
    updateSelectizeInput(session, "ontology_choice", choices=ontology_df$name, server=TRUE)
  })

  output$selected_var <- renderUI({
    sel <- input$table_rows_selected
    if (length(sel)) HTML(paste0("<strong>Selected:</strong> ", rv$data$variable[sel]))
    else HTML("<em>No row selected</em>")
  })

  # output$table <- renderDT({
  #   req(rv$data)
  #
  #   # Identify index of 'category' column
  #   category_col_index <- which(colnames(rv$data) == "category")
  #
  #   # All other columns should be disabled for editing
  #   disabled_cols <- setdiff(seq_along(rv$data), category_col_index) - 1
  #
  #
  #   datatable(
  #     rv$data,
  #     selection = "multiple",
  #     rownames = FALSE,
  #     editable = list(target = "cell", disable = list(columns = disabled_cols)),
  #     options = list(
  #       pageLength = 2000,
  #       scrollY = "400px",
  #       scrollX = TRUE,
  #       autoWidth = TRUE,
  #       stateSave = TRUE
  #     ),
  #     class = 'stripe nowrap'
  #   ) %>%
  #     formatStyle(
  #       'category_source',
  #       target = 'row',
  #       backgroundColor = styleEqual("manual", '#fff7e6')
  #     )
  # }, server = TRUE)

  output$table <- renderDT({
    req(rv$data)

    # Determine row status
    display_data <- rv$data
    display_data$row_status <- case_when(
      is.na(display_data$selected_ontology_id) ~ "Missing Annotation",
      is.na(display_data$category)             ~ "Missing Category",
      display_data$category_source == "manual" ~ "Manual",
      TRUE                                     ~ "OK"
    )

    datatable(
      display_data,
      selection = "multiple",
      rownames = FALSE,
      editable = list(target = "cell", disable = list(columns = which(names(display_data) != "category") - 1)),
      options = list(
        pageLength = 2000,
        scrollY = "400px",
        scrollX = TRUE,
        autoWidth = TRUE,
        stateSave = TRUE,
        columnDefs = list(list(visible = FALSE, targets = which(names(display_data) == "row_status") - 1))
      ),
      class = 'stripe nowrap'
    ) %>%
      formatStyle(
        'row_status',
        target = 'row',
        backgroundColor = styleEqual(
          c("Missing Annotation", "Missing Category", "Manual"),
          c("#fde0dc",             "#fff2cc",          "#fff7e6")
        )
      )
  }, server = TRUE)




  proxy <- dataTableProxy("table")

  observeEvent(input$apply_annotation, {
    req(rv$data, input$ontology_choice)
    rows <- input$table_rows_selected
    if (length(rows)) {
      sel <- ontology_df %>% filter(name == input$ontology_choice)
      if (nrow(sel)==1) {
        rv$data$selected_ontology_label[rows] <- sel$name
        rv$data$selected_ontology_id[rows]    <- sel$id
      }
    }
    replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
  })

  output$ontology_description <- renderUI({
    req(input$ontology_choice)
    term <- ontology_df %>% filter(name==input$ontology_choice)
    if (nrow(term)==1) {
      desc <- get_ols_description(term$id, strsplit(term$id,":")[[1]][1])
      if (!is.na(desc)) {
        div(style="border-left:4px solid #71255A;padding-left:10px;margin:10px 0;",
            strong("Description:"), br(), desc)
      }
    }
  })

  observeEvent(input$apply_categorization, {
    req(rv$data)
    if (all(is.na(rv$data$selected_ontology_id))) {
      showModal(modalDialog("Annotate at least one ID before categorizing.", easyClose=TRUE))
      return()
    }
    root_lvl <- as.numeric(input$cat_root_depth)
    rows <- input$table_rows_selected
    if (length(rows) == 0) {
      showModal(modalDialog("Please select one or more rows to categorize.", easyClose = TRUE))
      return()
    }

    # Only run categorization on selected rows
    data_subset <- rv$data[rows, , drop = FALSE]
    updated <- run_categorize_ontology(
      data          = data_subset,
      id_col        = "selected_ontology_id",
      ontology      = input$cat_ontology,
      root_level    = root_lvl,
      assign_label  = TRUE
    )

    # Only update selected rows
    rv$data$root_id[rows]         <- updated$root_id
    rv$data$root_label[rows]      <- updated$root_label
    rv$data$category[rows]        <- updated$category
    rv$data$category_source[rows] <- updated$category_source


    replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)

  })

  observeEvent(input$table_cell_edit, {
    info <- input$table_cell_edit
    print(info)  # debug

    i <- info$row
    j <- info$col
    v <- info$value

    colname <- names(rv$data)[j + 1]  # <-- FIX: adjust for zero-based index

    if (
      colname == "category" &&
      (is.na(rv$data$root_id[i]) || nchar(rv$data$root_id[i]) <= 1)
    ) {
      rv$data$category[i] <- v
      rv$data$category_source[i] <- "manual"
      replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
    } else {
      showModal(modalDialog(
        title = "Editing not allowed",
        "Only 'Unmapped' categories can be changed manually.",
        easyClose = TRUE
      ))
    }
  })




  output$save_csv <- downloadHandler(
    filename = function() paste0("annotated_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(rv$data, file, row.names=FALSE)
  )
}

# ---- Launch ----
shinyApp(ui, server)
