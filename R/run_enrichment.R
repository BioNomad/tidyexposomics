#' Run functional enrichment on omics features
#'
#' Runs over-representation analysis with \pkg{fenr} on features selected
#' from an `MultiAssayExperiment`-based exposomicset. Features are mapped to
#' gene identifiers, tested against a per-experiment universe (all features
#' in that experiment's `rowData`), p-values are adjusted across all tests,
#' and enriched terms are clustered by gene overlap (Jaccard similarity).
#'
#' @param exposomicset A `MultiAssayExperiment` (exposomicset) object with
#'   the upstream analysis results stored in its metadata.
#' @param feature_type Character. Which features to test:
#'   \describe{
#'     \item{`"degs"`}{Differentially abundant features passing the p-value
#'       and log fold-change thresholds.}
#'     \item{`"degs_robust"`}{Features with a stability score above
#'       `score_threshold` from the sensitivity analysis.}
#'     \item{`"omics"`}{Features supplied via `variable_map`.}
#'     \item{`"factor_features"`}{Features loading on multi-omics
#'       integration factors (see `factor_type`).}
#'     \item{`"degs_cor"`, `"omics_cor"`, `"factor_features_cor"`}{Features
#'       from the corresponding correlation results, tested per
#'       `category`.}
#'   }
#' @param score_col Character. Column in the feature stability table used
#'   for `"degs_robust"`. Default `"stability_score"`.
#' @param score_threshold Numeric or `NULL`. Minimum stability score. If
#'   `NULL`, the threshold stored in the sensitivity analysis is used.
#' @param variable_map A data frame with a `variable` column (and `exp_name`,
#'   `category` where applicable) mapping variables to features. Required
#'   when `feature_type` is `"omics"`.
#' @param factor_type Character. Which factor feature set to use for
#'   `"factor_features"`: `"common_top_factor_features"` or
#'   `"top_factor_features"`.
#' @param feature_col Character. Column in the feature metadata
#'   (`rowData`) holding the identifier used for enrichment, e.g. gene
#'   symbols. Default `"feature"`.
#' @param deg_pval_col Character. Column with adjusted p-values in the
#'   differential abundance table. Default `"adj.P.Val"`.
#' @param deg_pval_threshold Numeric. P-value cutoff for calling DEGs
#'   (strictly less than). Default `0.05`.
#' @param deg_logfc_col Character. Column with log fold-changes in the
#'   differential abundance table. Default `"logFC"`.
#' @param deg_logfc_threshold Numeric. Absolute log fold-change cutoff for
#'   calling DEGs (strictly greater than). Default `log2(1.5)`.
#' @param db Character. Annotation database: `"GO"`, `"KEGG"`, or
#'   `"Reactome"`.
#' @param species Character or `NULL`. Species designation passed to
#'   \pkg{fenr}. Required for `db = "GO"`, and used for Reactome.
#' @param fenr_col Character. Name of the feature column in the
#'   \pkg{fenr} term mapping to match against. Default `"gene_symbol"`.
#' @param padj_method Character. Method passed to [stats::p.adjust()].
#'   Default `"fdr"`.
#' @param pval_thresh Numeric. Adjusted p-value cutoff for retained terms.
#'   Default `0.1`.
#' @param min_set,max_set Numeric. Terms are kept only if the number of
#'   selected features annotated to them is greater than `min_set` and less
#'   than `max_set`. Defaults `3` and `800`.
#' @param clustering_approach Character. Clustering method used to group
#'   enriched terms from the Jaccard similarity matrix. Default `"diana"`.
#' @param action Character. `"add"` stores the results in the object's
#'   metadata and returns the object; `"get"` returns the enrichment
#'   results directly.
#'
#' @return If `action = "get"`, a tibble of enrichment results with
#'   `exp_name`, adjusted p-values (`padj`), a `go_group` cluster
#'   assignment and, for `_cor` feature types, `category`. If
#'   `action = "add"`, the input `exposomicset` with results stored in
#'   `metadata(exposomicset)$enrichment[[feature_type]]` and the run
#'   parameters recorded in `metadata(exposomicset)$summary$steps`.
#'
#' @details
#' Term data are downloaded only when there are features to test. For
#' `db = "GO"`, \pkg{fenr}'s default annotation URL is overridden and
#' `fenr::fetch_go()` is bypassed, because both rely on endpoints that are
#' no longer available. Fewer than two selected genes in an
#' experiment/category skips that test; fewer than two enriched terms skips
#' clustering and sets `go_group` to `NA`.
#'
#' @examples
#' \dontrun{
#' # Enrichment of DEGs against GO (biological process terms)
#' mae <- run_enrichment(
#'     exposomicset = mae,
#'     feature_type = "degs",
#'     db = "GO",
#'     species = "hsapiens"
#' )
#'
#' # Return results without modifying the object
#' res <- run_enrichment(
#'     exposomicset = mae,
#'     feature_type = "degs_robust",
#'     db = "KEGG",
#'     action = "get"
#' )
#' }
#'
#' @importFrom MultiAssayExperiment metadata experiments
#' @importFrom SummarizedExperiment rowData
#' @importFrom rlang .data
#' @importFrom stats p.adjust as.dist
#' @importFrom Matrix sparseMatrix tcrossprod rowSums
#' @import fenr
#' @export
run_enrichment <- function(
  exposomicset,
  feature_type = c(
      "degs",
      "degs_robust",
      "omics",
      "factor_features",
      "degs_cor",
      "omics_cor",
      "factor_features_cor"
  ),
  score_col = "stability_score",
  score_threshold = NULL,
  variable_map = NULL,
  factor_type = c(
      "common_top_factor_features",
      "top_factor_features"
  ),
  feature_col = "feature",
  deg_pval_col = "adj.P.Val",
  deg_pval_threshold = 0.05,
  deg_logfc_col = "logFC",
  deg_logfc_threshold = log2(1.5),
  db = c("GO", "KEGG", "Reactome"),
  species = NULL,
  fenr_col = "gene_symbol",
  padj_method = "fdr",
  pval_thresh = 0.1,
  min_set = 3,
  max_set = 800,
  clustering_approach = "diana",
  action = c("add", "get")
) {
    db <- match.arg(db)
    factor_type <- match.arg(factor_type)
    feature_type <- match.arg(feature_type)
    action <- match.arg(action)

    if (db == "GO" && is.null(species)) {
        stop("Please specify a species designation for GO.")
    }

    enr_res <- .select_features(
        exposomicset = exposomicset,
        feature_type = feature_type,
        fdata = if (feature_type != "degs") {
            .get_fdata(exposomicset, feature_col)
        },
        feature_col = feature_col,
        deg_pval_col = deg_pval_col,
        deg_pval_threshold = deg_pval_threshold,
        deg_logfc_col = deg_logfc_col,
        deg_logfc_threshold = deg_logfc_threshold,
        score_col = score_col,
        score_threshold = score_threshold,
        factor_type = factor_type,
        variable_map = variable_map
    ) |>
        .run_enr_core(
            exposomicset = exposomicset,
            # lazy: only downloaded if there are genes to test
            term_data = .fetch_term_data(db, species),
            feature_col = feature_col,
            fenr_col = fenr_col,
            padj_method = padj_method,
            pval_thresh = pval_thresh,
            min_set = min_set,
            max_set = max_set
        ) |>
        .group_enr_res(clustering_approach = clustering_approach)

    if (action == "get") {
        return(enr_res)
    }

    # --- Store results ---
    all_metadata <- MultiAssayExperiment::metadata(exposomicset)
    all_metadata$enrichment[[feature_type]] <- enr_res
    MultiAssayExperiment::metadata(exposomicset) <- all_metadata

    MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
        MultiAssayExperiment::metadata(exposomicset)$summary$steps,
        list(
            run_enrichment = list(
                timestamp = Sys.time(),
                params = list(
                    feature_type = feature_type,
                    score_col = score_col,
                    score_threshold = score_threshold,
                    variable_map = variable_map,
                    factor_type = factor_type,
                    feature_col = feature_col,
                    deg_pval_col = deg_pval_col,
                    deg_pval_threshold = deg_pval_threshold,
                    deg_logfc_col = deg_logfc_col,
                    deg_logfc_threshold = deg_logfc_threshold,
                    db = db,
                    species = species,
                    fenr_col = fenr_col,
                    padj_method = padj_method,
                    pval_thresh = pval_thresh,
                    clustering_approach = clustering_approach
                ),
                notes = paste0(
                    "Performed ", db, " enrichment on ", feature_type,
                    " features."
                )
            )
        )
    )

    exposomicset
}

# --- Fetch term data ------------
#' @keywords internal
#' @noRd
.fetch_term_data <- function(db, species) {
    if (db == "GO") {
        # fenr's default GO URL is dead
        options(
            GO_ANNOTATION_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN"
        )
    }

    tryCatch(
        switch(db,
            GO = list(
                # bypass fenr::fetch_go(), whose assert_species() hits a dead URL
                terms = fenr:::fetch_go_terms(
                    use_cache = TRUE,
                    on_error = "stop"
                ),
                mapping = fenr:::fetch_go_genes_go(
                    species = species,
                    use_cache = TRUE,
                    on_error = "stop"
                )
            ),
            KEGG = fenr::fetch_kegg(),
            Reactome = fenr::fetch_reactome(
                species = species,
                use_cache = TRUE
            )
        ),
        error = \(e) {
            message(e$message)
            message(db, " term data failed to download")
            NULL
        }
    )
}

# --- Prepare term data for one universe  ------------
#' @keywords internal
#' @noRd
.prep_fenr <- function(term_data, universe_genes, fenr_col) {
    tryCatch(
        fenr::prepare_for_enrichment(
            terms = term_data$terms,
            mapping = term_data$mapping,
            all_features = universe_genes,
            feature_name = fenr_col
        ),
        error = \(e) {
            message(e$message)
            NULL
        }
    )
}

# --- Run enrichment for one gene set ------------
#' @keywords internal
#' @noRd
.run_fenr <- function(selected_genes, universe_genes, terms_obj) {
    tryCatch(
        fenr::functional_enrichment(
            feat_all = universe_genes,
            feat_sel = selected_genes,
            term_data = terms_obj
        ) |>
            tibble::as_tibble(),
        error = \(e) {
            message(e$message)
            tibble::tibble()
        }
    )
}

# --- Minimal fdata: only the columns needed for joins ------------
#' @keywords internal
#' @noRd
.get_fdata <- function(exposomicset, feature_col) {
    exposomicset |>
        pivot_feature() |>
        dplyr::select(
            dplyr::any_of(c(".feature", ".exp_name", feature_col))
        ) |>
        dplyr::mutate(
            dplyr::across(dplyr::all_of(feature_col), as.character)
        ) |>
        dplyr::distinct()
}

# --- Map feature IDs to gene IDs ------------
#' @keywords internal
#' @noRd
.join_fdata <- function(x, fdata, feature_col) {
    x |>
        dplyr::select(dplyr::any_of(c("feature", "exp_name", "category"))) |>
        dplyr::inner_join(
            fdata,
            by = c("feature" = ".feature", "exp_name" = ".exp_name")
        ) |>
        dplyr::select(
            dplyr::any_of(c("exp_name", "category")),
            gene = dplyr::all_of(feature_col)
        ) |>
        dplyr::distinct()
}

# --- Plain character vector, no NA/empty ------------
#' @keywords internal
#' @noRd
.as_chr <- function(x) {
    x |>
        unlist(use.names = FALSE) |>
        as.character() |>
        (\(v) unique(v[!is.na(v) & nzchar(v)]))()
}

# --- Selected genes for every experiment (and category) at once ------------
#' @return tibble with `exp_name`, `gene` and, for `_cor` types, `category`
#' @keywords internal
#' @noRd
.select_features <- function(
  exposomicset,
  feature_type,
  fdata,
  feature_col,
  deg_pval_col,
  deg_pval_threshold,
  deg_logfc_col,
  deg_logfc_threshold,
  score_col,
  score_threshold,
  factor_type,
  variable_map
) {
    md <- exposomicset@metadata

    switch(feature_type,
        "degs" = md |>
            purrr::pluck("differential_analysis", "differential_abundance") |>
            dplyr::filter(
                .data[[deg_pval_col]] < deg_pval_threshold,
                abs(.data[[deg_logfc_col]]) > deg_logfc_threshold
            ) |>
            dplyr::select(exp_name, gene = dplyr::all_of(feature_col)) |>
            dplyr::mutate(gene = as.character(gene)) |>
            dplyr::distinct(),
        "degs_robust" = {
            sens <- purrr::pluck(
                md, "differential_analysis", "sensitivity_analysis"
            )

            sens |>
                purrr::pluck("feature_stability") |>
                dplyr::filter(
                    .data[[score_col]] > (score_threshold %||% sens$score_thresh)
                ) |>
                .join_fdata(fdata, feature_col)
        },
        "factor_features" = md |>
            purrr::pluck(
                "differential_analysis", "multiomics_integration", factor_type
            ) |>
            .join_fdata(fdata, feature_col),
        "omics" = variable_map |>
            dplyr::rename(feature = variable) |>
            .join_fdata(fdata, feature_col),
        "degs_cor" = ,
        "omics_cor" = ,
        "factor_features_cor" = md |>
            purrr::pluck("correlation", sub("_cor$", "", feature_type)) |>
            .join_fdata(fdata, feature_col)
    )
}

# --- Run fenr per experiment (and category), adjust, filter ------------
#' @keywords internal
#' @noRd
.run_enr_core <- function(
  sel,
  exposomicset,
  term_data,
  feature_col,
  fenr_col,
  padj_method,
  pval_thresh,
  min_set,
  max_set
) {
    if (nrow(sel) == 0L || is.null(term_data)) {
        return(tibble::tibble())
    }

    has_cat <- "category" %in% names(sel)

    enr_res <- sel |>
        split(sel$exp_name) |>
        purrr::imap(\(d_exp, exp) {
            universe <- SummarizedExperiment::rowData(
                exposomicset[[exp]]
            )[[feature_col]] |>
                .as_chr()

            # one prepare_for_enrichment() per experiment, reused across categories
            terms_obj <- if (length(universe) > 0L) {
                .prep_fenr(term_data, universe, fenr_col)
            }

            if (is.null(terms_obj)) {
                return(NULL)
            }

            (if (has_cat) split(d_exp, d_exp$category) else list(d_exp)) |>
                purrr::map(\(chunk) {
                    genes <- unique(chunk$gene)

                    if (length(genes) < 2L) {
                        return(NULL)
                    }

                    .run_fenr(genes, universe, terms_obj) |>
                        dplyr::mutate(
                            exp_name = exp,
                            category = if (has_cat) chunk$category[[1]]
                        )
                })
        }) |>
        purrr::list_flatten() |>
        dplyr::bind_rows()

    if (nrow(enr_res) == 0L) {
        return(tibble::tibble())
    }

    # Multiple-hypothesis correction after all tests are done
    enr_res |>
        dplyr::mutate(padj = stats::p.adjust(p_value, method = padj_method)) |>
        dplyr::filter(
            padj < pval_thresh,
            n_with_sel > min_set,
            n_with_sel < max_set
        )
}

# --- Pairwise Jaccard similarity via one sparse matrix product ------------
#' @param gene_list Named list of character vectors (one per term)
#' @return Square numeric similarity matrix with term names as dimnames
#' @keywords internal
#' @noRd
.jaccard_mat <- function(gene_list) {
    genes <- unlist(gene_list, use.names = FALSE)

    m <- Matrix::sparseMatrix(
        i = rep(seq_along(gene_list), lengths(gene_list)),
        j = match(genes, genes),
        x = 1,
        dims = c(length(gene_list), length(genes)),
        dimnames = list(names(gene_list), NULL)
    )

    inter <- as.matrix(Matrix::tcrossprod(m))
    n <- Matrix::rowSums(m)

    inter / (outer(n, n, "+") - inter)
}

# --- Group Enrichment Terms ------------
#' @keywords internal
#' @noRd
.group_enr_res <- function(enr_res, clustering_approach = "diana") {
    if (nrow(enr_res) == 0L ||
        !all(c("term_name", "ids") %in% names(enr_res))) {
        return(enr_res)
    }

    if (nrow(enr_res) < 2L) {
        message("GO clustering skipped - fewer than 2 enriched terms.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    message("Determining Number of GO Term Clusters...")

    gene_list <- enr_res |>
        dplyr::distinct(term_name, ids) |>
        (\(d) strsplit(d$ids, ",") |>
            lapply(trimws) |>
            split(d$term_name))() |>
        lapply(\(x) unique(unlist(x)))

    if (length(gene_list) < 2L) {
        message("GO clustering skipped - fewer than 2 unique terms.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    go_groups <- gene_list |>
        .jaccard_mat() |>
        .cluster_mat(clustering_approach = clustering_approach) |>
        (\(x) data.frame(
            term_name = names(x),
            go_group = paste0("Group ", as.numeric(x))
        ))()

    if (is.null(go_groups) || nrow(go_groups) == 0L) {
        message("GO clustering failed - no clusters found.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    dplyr::inner_join(enr_res, go_groups, by = "term_name")
}
