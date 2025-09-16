#' Internal - load ontologies
#'
#' @keywords internal
#' @param use_demo Logical; if TRUE, load packaged demo objects \code{hpo}, \code{ecto}, \code{chebi}.
#' @param ... Optional named overrides: \code{hpo}, \code{ecto}, \code{chebi}.
#' @return A list with elements \code{hpo}, \code{ecto}, \code{chebi}.
#' @importFrom utils data
.load_ontologies <- function(use_demo = TRUE, ...) {
    dots <- list(...)
    has <- function(nm) !is.null(dots[[nm]])

    if (use_demo) {
        if (exists("load_annotation_data", mode = "function")) {
            tidyexposomics::load_annotation_data()
        } else {
            utils::data(list = c("hpo", "ecto", "chebi"), package = "tidyexposomics", envir = environment())
        }
    }

    list(
        hpo   = if (has("hpo")) dots$hpo else get("hpo", inherits = TRUE),
        ecto  = if (has("ecto")) dots$ecto else get("ecto", inherits = TRUE),
        chebi = if (has("chebi")) dots$chebi else get("chebi", inherits = TRUE)
    )
}

#' Internal - categorizes to ontology roots / depth
#'
#' @keywords internal
#' @param data Data frame with an ID column.
#' @param id_col Column name containing ontology term IDs.
#' @param ontologyDF Processed ontology data.frame (with list columns for relationships and \code{depth}).
#' @param root_level Either a character vector of explicit root IDs, a numeric depth, or (default) top-level roots.
#' @param assign_label If TRUE, return input with new columns; otherwise return assigned IDs.
#' @return Data frame with \code{root_id}, \code{root_label}, \code{category}, \code{category_source} (if \code{assign_label}).
#' @importFrom dplyr mutate select any_of
#' @importFrom stats setNames
.run_categorize_ontology <- function(data, id_col, ontologyDF, root_level = 0, assign_label = TRUE) {
    id_to_ancestors <- stats::setNames(ontologyDF$ancestors, ontologyDF$id)
    id_to_depth <- stats::setNames(ontologyDF$depth, ontologyDF$id)
    id_to_name <- stats::setNames(ontologyDF$name, ontologyDF$id)

    assign_to_root <- switch(class(root_level),
        "character" = {
            root_nodes <- root_level
            function(term_id) {
                anc <- c(term_id, id_to_ancestors[[term_id]])
                matched <- intersect(anc, root_nodes)
                if (length(matched)) matched[1] else NA_character_
            }
        },
        "numeric" = {
            function(term_id) {
                anc <- c(term_id, id_to_ancestors[[term_id]])
                depths <- id_to_depth[anc]
                idxs <- which(depths == root_level)
                if (length(idxs)) anc[idxs[1]] else NA_character_
            }
        },
        {
            root_nodes <- ontologyDF$id[lengths(ontologyDF$parents) == 0]
            function(term_id) {
                anc <- c(term_id, id_to_ancestors[[term_id]])
                matched <- intersect(anc, root_nodes)
                if (length(matched)) matched[1] else NA_character_
            }
        }
    )

    term_ids <- data[[id_col]]
    assigned_ids <- vapply(term_ids, assign_to_root, character(1))
    assigned_labels <- id_to_name[assigned_ids]

    if (assign_label) {
        dplyr::select(data, -dplyr::any_of(c("root_id", "root_label", "category", "category_source"))) |>
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

#' Internal - fetch OLS description for a term
#'
#' @keywords internal
#' @param ontology_id Term ID (e.g., \code{"HP:0004322"}).
#' @param ontology_prefix Prefix (e.g., \code{"HP"}).
#' @return A character scalar description or \code{NA_character_}.
#' @importFrom utils URLencode
#' @importFrom httr GET timeout status_code
#' @importFrom jsonlite fromJSON
.get_ols_description <- function(ontology_id, ontology_prefix) {
    iri <- paste0("http://purl.obolibrary.org/obo/", gsub(":", "_", ontology_id))
    double_encode <- function(url) {
        tmp <- utils::URLencode(url, reserved = TRUE)
        gsub("%", "%25", tmp, fixed = TRUE)
    }
    enc <- double_encode(iri)
    url <- paste0(
        "https://www.ebi.ac.uk/ols4/api/ontologies/",
        tolower(ontology_prefix), "/terms/", enc
    )

    if (!.can_use_network()) {
        return(NA_character_)
    }

    res <- try(httr::GET(url, httr::timeout(3)), silent = TRUE)
    if (inherits(res, "try-error") || httr::status_code(res) != 200) {
        return(NA_character_)
    }
    j <- jsonlite::fromJSON(rawToChar(res$content))
    desc <- j$description
    if (!is.null(desc) && length(desc) > 0 && nzchar(desc[[1]])) {
        desc[[1]]
    } else if (!is.null(j$annotation$description) && length(j$annotation$description) > 0) {
        j$annotation$description[[1]]
    } else {
        NA_character_
    }
}

#' Internal - should we allow network access now?
#'
#' Controlled by the R option \code{tidyexposomics.allow_network} and \code{interactive()}.
#'
#' @keywords internal
#' @return Logical scalar.
.can_use_network <- function() {
    getOption("tidyexposomics.allow_network", TRUE) && interactive()
}
