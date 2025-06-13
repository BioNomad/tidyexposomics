#' Categorize ontology terms by root-level ancestor (Optimized)
#'
#' This function maps ontology term IDs to higher-level ancestor categories by
#' traversing the hierarchy of terms in the Chemical Entities of
#' Biological Interest (ChEBI), Human Phenotype Ontology (HPO),
#' and Experimental Conditions Ontology (ECTO).
#'
#' @param data A data frame containing ontology term IDs to categorize.
#' @param id_col A string specifying the name of the column in `data` that contains ontology term IDs.
#' @param root_level Either `"top"` (default) to assign each term to its topmost ancestor,
#'   a numeric depth to assign to the closest ancestor at or below that depth, or
#'   a character vector of term IDs to use as custom roots.
#' @param assign_label Logical; if `TRUE` (default), returns labels for both original and root terms.
#'
#' @return If `assign_label = TRUE`, returns a data frame with additional columns:
#'   \describe{
#'     \item{root_id}{ID of the assigned root ancestor}
#'     \item{root_label}{Label of the assigned root ancestor}
#'     \item{category}{Alias of `root_label`, convenient for downstream use}
#'   }
#'   Otherwise, returns a character vector of root IDs.
#'
#' @export
run_categorize_ontology <- function(
    data,
    id_col,
    ontology,
    root_level = "top",
    assign_label = TRUE
) {

  message("Downloading ", ontology, " ontology data.")
  ontology_df <- switch(
    ontology,
    "hpo" = readRDS("../data/hpo.rds"),
    "ecto" = readRDS("../data/ecto.rds"),
    "chebi" = readRDS("../data/chebi.rds"),
    stop("Invalid ontology specified. Choose from 'hpo', 'ecto', or 'chebi'.")
  )
  # ontology_df <- switch(
  #   ontology,
  #   "hpo" = ontologyIndex::get_ontology(
  #       file = url("https://purl.obolibrary.org/obo/hp.obo")) |>
  #     as.data.frame(),
  #   "ecto" = ontologyIndex::get_ontology(
  #       file = url("https://purl.obolibrary.org/obo/ecto.obo")) |>
  #     as.data.frame(),
  #   "chebi" = ontologyIndex::get_ontology(
  #     file = url("https://purl.obolibrary.org/obo/chebi.obo")) |>
  #     as.data.frame(),
  #   stop("Invalid ontology specified. Choose from 'hpo', 'ecto', or 'chebi'.")
  # )

  # Helper to fix relationship columns if needed
  fix_rel_cols <- function(x) {
    if (!is.list(x)) strsplit(x, ";\\s*") else x
  }

  ontology_df <- ontology_df |>
    dplyr::mutate(
      ancestors = fix_rel_cols(ancestors),
      parents = fix_rel_cols(parents),
      children = fix_rel_cols(children),
      depth = lengths(ancestors)
    )

  # Precompute lookup tables for speed
  id_to_ancestors <- setNames(ontology_df$ancestors, ontology_df$id)
  id_to_depth <- setNames(ontology_df$depth, ontology_df$id)
  id_to_name <- setNames(ontology_df$name, ontology_df$id)

  # Define root assignment logic
  assign_to_root <- switch(
    class(root_level),
    character = {
      root_nodes <- root_level
      function(term_id) {
        anc <- c(term_id, id_to_ancestors[[term_id]])
        matched <- intersect(anc, root_nodes)
        if (length(matched)) matched[1] else NA_character_
      }
    },
    numeric = {
      function(term_id) {
        anc <- c(term_id, id_to_ancestors[[term_id]])
        depths <- id_to_depth[anc]
        depths <- depths[!is.na(depths)]
        eligible <- depths[depths <= root_level]
        if (length(eligible)) names(eligible)[which.max(eligible)] else NA_character_
      }
    },
    default = {
      root_nodes <- ontology_df$id[lengths(ontology_df$parents) == 0]
      function(term_id) {
        anc <- c(term_id, id_to_ancestors[[term_id]])
        matched <- intersect(anc, root_nodes)
        if (length(matched)) matched[1] else NA_character_
      }
    }
  )

  # Apply categorization
  term_ids <- data[[id_col]]
  assigned_root_ids <- vapply(term_ids, assign_to_root, character(1))

  if (assign_label) {
    dplyr::bind_cols(
      data,
      tibble::tibble(
        root_id = assigned_root_ids,
        root_label = unname(id_to_name[assigned_root_ids])
      )
    ) |>
      dplyr::mutate(category = root_label)
  } else {
    assigned_root_ids
  }
}

run_categorize_mult_ontology <- function(
    data,
    pattern,
    root_level,
    ontology,
    id_col = "selected_ontology_id"
) {
  # Input validation
  if (!(length(pattern) == length(root_level) && length(root_level) == length(ontology))) {
    stop("`pattern`, `root_level`, and `ontology` must be the same length.")
  }

  # Zip the arguments and apply in parallel
  purrr::pmap_dfr(
    list(pattern = pattern, root_level = root_level, ontology = ontology),
    function(pattern, root_level, ontology) {
      data |>
        dplyr::filter(grepl(pattern, .data[[id_col]])) |>
        run_categorize_ontology(
          id_col = id_col,
          root_level = root_level,
          ontology = ontology,
          assign_label = TRUE
        )
    }
  )
}



