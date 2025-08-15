#' Perform enrichment analysis on selected features from a expomicset object
#'
#' This function performs enrichment analysis using selected features derived
#' from differential expression, correlation analysis,
#' or multi-omics factor features across experiments in an `expomicset`.
#' It supports multiple enrichment databases (e.g., GO, KEGG, Reactome),
#'  applies FDR correction, and optionally clusters GO terms by Jaccard overlap.
#'
#' @param expomicset An `expomicset` (a `MultiAssayExperiment`
#' object with metadata) containing omics and metadata.
#' @param feature_type Character string indicating the feature source.
#' One of `"degs"`, `"degs_robust"`,
#'   `"omics"`, `"factor_features"`, `"degs_cor"`, `"omics_cor"`,
#'    or `"factor_features_cor"`.
#' @param score_col Column name used for stability score filtering
#' (only for `degs_robust`).
#' @param score_threshold Optional numeric threshold for filtering
#' stability scores. If `NULL`, the default
#'   threshold stored in the metadata will be used.
#' @param variable_map A data frame with `exp_name` and `feature` columns,
#'  used when `feature_type = "omics"`.
#' @param factor_type Character string for selecting factor features:
#'  `"common_top_factor_features"` or
#'   `"top_factor_features"`.
#' @param feature_col The name of the feature column used to extract
#'  gene identifiers.
#' @param deg_pval_col Column name for adjusted p-values from DEG analysis.
#' @param deg_pval_threshold Threshold to select significant DEGs
#' (default: 0.05).
#' @param deg_logfc_col Column name for log-fold changes from DEG analysis.
#' @param deg_logfc_threshold Threshold to select DEGs by absolute logFC
#' (default: `log2(1.5)`).
#' @param db Enrichment database to use. One of `"GO"`,
#'  `"KEGG"`, `"Reactome"`.
#' @param species Species name (required for GO enrichment,
#' e.g., `"Homo sapiens"`). Ignored for other databases.
#' @param fenr_col Column name for gene IDs used by `fenr`
#'  (e.g., `"gene_symbol"`).
#' @param padj_method Method for p-value adjustment (default: `"fdr"`).
#' @param pval_thresh Adjusted p-value threshold for filtering enriched terms
#' (default: 0.1).
#' @param min_set Minimum number of selected genes overlapping an enriched term
#' (default: 3).
#' @param max_set Maximum number of selected genes overlapping an enriched term
#' (default: 800).
#' @param clustering_approach Clustering method for GO term grouping.
#'  Defaults to `"diana"`.
#' @param action Either `"add"` to store results in the object's metadata
#'  or `"get"` to return results as a data frame.
#'
#' @return If `action = "add"`, returns the modified `ExpOmicSet` with
#' enrichment results added to metadata.
#' If `action = "get"`, returns a `data.frame` of enrichment results
#' with GO term clusters (if applicable).
#'
#' @details
#' The function identifies selected features based on the chosen
#' `feature_type`, determines the gene universe
#' for each experiment, and performs enrichment analysis using the
#' `fenr` package. Results are adjusted for
#' multiple testing and optionally clustered by gene set overlap (for GO terms).
#'
#' If `feature_type` includes correlation-based results
#' (ending in `_cor`), enrichment is performed for each
#' exposure category separately.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 30,
#'     return_mae = TRUE
#' )
#'
#' # perform differential abundance analysis
#' mae <- run_differential_abundance(
#'     expomicset = mae,
#'     formula = ~ smoker + sex,
#'     abundance_col = "counts",
#'     method = "limma_voom",
#'     action = "add"
#' )
#'
#' # perform enrichment analysis
#' mae <- run_enrichment(
#'     expomicset = mae,
#'     feature_type = "degs",
#'     feature_col = "symbol",
#'     species = "goa_human",
#'     deg_logfc_threshold = log2(1),
#'     deg_pval_col = "P.Value",
#'     deg_pval_threshold = 0.2,
#'     action = "add"
#' )
#'
#' @importFrom MultiAssayExperiment metadata experiments
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr filter inner_join pull mutate summarise bind_rows
#' select distinct arrange
#' @importFrom tidyr unnest pivot_wider
#' @importFrom purrr pluck map map2
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_split str_trim
#' @importFrom rlang sym set_names .data
#' @importFrom stats p.adjust
#' @import fenr
#' @export
run_enrichment <- function(
    expomicset,
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
    action = "add") {
    # require(fenr)
    db <- match.arg(db)
    factor_type <- match.arg(factor_type)
    feature_type <- match.arg(feature_type)

    # Grab all feature meta data
    fdata <- expomicset |>
        pivot_feature()

    if (feature_type %in% c(
        "degs",
        "degs_robust",
        "omics",
        "factor_features"
    )) {
        enr_res <- .run_feature_enr(
            expomicset = expomicset,
            feature_type = feature_type,
            db = db,
            fdata = fdata,
            species = species,
            fenr_col = fenr_col,
            feature_col = feature_col,
            deg_pval_col = deg_pval_col,
            deg_pval_threshold = deg_pval_threshold,
            deg_logfc_col = deg_logfc_col,
            deg_logfc_threshold = deg_logfc_threshold,
            score_col = score_col,
            score_threshold = score_threshold,
            factor_type = factor_type,
            variable_map = variable_map,
            padj_method = padj_method,
            pval_thresh = pval_thresh,
            min_set = min_set,
            max_set = max_set
        )
    }

    if (feature_type %in% c(
        "degs_cor",
        "omics_cor",
        "factor_features_cor"
    )) {
        enr_res <- .run_feature_cor_enr(
            expomicset = expomicset,
            feature_type = feature_type,
            db = db,
            fdata = fdata,
            species = species,
            fenr_col = fenr_col,
            feature_col = feature_col,
            deg_pval_col = deg_pval_col,
            deg_pval_threshold = deg_pval_threshold,
            deg_logfc_col = deg_logfc_col,
            deg_logfc_threshold = deg_logfc_threshold,
            score_col = score_col,
            score_threshold = score_threshold,
            factor_type = factor_type,
            variable_map = variable_map,
            padj_method = padj_method,
            pval_thresh = pval_thresh,
            min_set = min_set,
            max_set = max_set
        )
    }

    # --- Group Terms ------------
    enr_res <- .group_enr_res(
        enr_res,
        clustering_approach = clustering_approach
    )

    # --- Store or return ---
    if (action == "add") {
        all_metadata <- MultiAssayExperiment::metadata(expomicset)
        all_metadata$enrichment[[feature_type]] <- enr_res
        MultiAssayExperiment::metadata(expomicset) <- all_metadata

        step_record <- list(
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
                    "Performed ",
                    db,
                    " enrichment on ",
                    feature_type,
                    " features."
                )
            )
        )

        MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(expomicset)$summary$steps,
            step_record
        )

        return(expomicset)
    } else if (action == "get") {
        return(enr_res)
    } else {
        stop("Invalid action. Choose from 'add' or 'get'.")
    }
}

# --- Run Fenr Enrichment Function ---------
#' @noRd
#'
#' @title Internal wrapper to run functional enrichment using `fenr`
#'
#' @description
#' Helper function to fetch database-specific term data,
#' prepare term mappings,
#' and perform enrichment analysis using the `fenr` package.
#'
#' @param selected_genes Character vector of selected gene identifiers.
#' @param universe_genes Character vector of all background/universe
#' gene identifiers.
#' @param db Character string specifying the database.
#' One of `"GO"`, `"KEGG"`, `"Reactome"`.
#' @param species Optional character string specifying species name
#' (required for GO).
#' @param fenr_col Column name in the term mappings corresponding to gene IDs
#' (default `"gene_symbol"`).
#'
#' @return A data frame of enrichment results from
#' `fenr::functional_enrichment()`, or `NULL` if enrichment fails.
#'
#' @importFrom MultiAssayExperiment metadata experiments
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr filter inner_join pull mutate summarise
#' bind_rows select distinct arrange
#' @importFrom tidyr unnest pivot_wider
#' @importFrom purrr pluck map map2
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_split str_trim
#' @importFrom rlang sym set_names .data
#' @importFrom stats p.adjust
#' @import fenr
#' @keywords internal
.run_fenr <- function(
    selected_genes,
    universe_genes,
    db = c("GO", "KEGG", "Reactome"),
    species = NULL,
    fenr_col = "gene_symbol") {
    # require(fenr)

    db <- match.arg(db)

    # Fetch functional terms + mapping
    fetch_fun <- switch(db,
        GO = fenr::fetch_go,
        KEGG = fenr::fetch_kegg,
        Reactome = fenr::fetch_reactome
    )

    # Handle species if required
    if (db == "GO" && is.null(species)) {
        stop("Please specify a species designation for GO.")
    }

    if (db == "GO") {
        term_data <- fetch_fun(species = species)
    } else {
        term_data <- fetch_fun()
    }

    # Prepare for enrichment
    # terms_obj <- fenr::prepare_for_enrichment(
    #   terms = term_data$terms,
    #   mapping = term_data$mapping,
    #   all_features = universe_genes,
    #   feature_name = fenr_col
    # )

    # Prepare for enrichment
    terms_obj <- tryCatch(
        {
            # may fail if there are no mappings
            fenr::prepare_for_enrichment(
                terms = term_data$terms,
                mapping = term_data$mapping,
                all_features = universe_genes,
                feature_name = fenr_col
            )
        },
        error = function(e) {
            message(e$message)
            return(NULL)
        }
    )


    # Run enrichment

    if (!is.null(terms_obj)) {
        enrichment_results <- fenr::functional_enrichment(
            feat_all = universe_genes,
            feat_sel = selected_genes,
            term_data = terms_obj
        )
    } else {
        enrichment_results <- NULL
    }

    return(enrichment_results)
}

# --- Run Enrichment on Features ------------
#' Run enrichment for feature-based selections across experiments
#'
#' This internal helper performs enrichment analysis for a given
#' \code{feature_type}
#' across all assays/experiments in an \code{expomicset}
#' (a \code{MultiAssayExperiment}).
#' It selects features per experiment according to \code{feature_type},
#' computes the gene universe from \code{rowData()}, runs \code{.run_fenr()},
#' then applies multiple testing correction and size filtering.
#'
#' Supported \code{feature_type} values (no \code{*_cor} here):
#' \itemize{
#'   \item \code{"degs"} - select DEGs by adjusted p-value and absolute value
#'   of the logFC
#'   \item \code{"degs_robust"} - select features using
#'   sensitivity/stability scores
#'   \item \code{"factor_features"} - select features from multi-omics
#'    integration
#'   \item \code{"omics"} - select features using a \code{variable_map}
#' }
#'
#' @param expomicset A \code{MultiAssayExperiment} with experiments and
#'  metadata.
#' @param feature_type Character; one of \code{"degs"}, \code{"degs_robust"},
#'   \code{"factor_features"}, \code{"omics"}.
#' @param db Enrichment database: one of \code{"GO"}, \code{"KEGG"},
#' \code{"Reactome"}.
#' @param species Species identifier required by the underlying database
#' (used by \code{.run_fenr}).
#' @param fenr_col Column name of gene IDs expected by the enrichment backend
#'  (e.g., \code{"gene_symbol"}).
#' @param feature_col Column in \code{rowData()} providing gene identifiers
#' for the universe and mapping.
#' @param deg_pval_col Column of adjusted p-values in DEG results
#'  (for \code{feature_type = "degs"}).
#' @param deg_pval_threshold Numeric threshold for adjusted p-values in
#' DEG selection.
#' @param deg_logfc_col Column of log-fold change values in DEG results
#' (for \code{feature_type = "degs"}).
#' @param deg_logfc_threshold Numeric threshold for absolute log-fold change
#' in DEG selection.
#' @param score_col Column name for stability/sensitivity score
#'  (for \code{"degs_robust"}).
#' @param score_threshold Optional numeric threshold; if \code{NULL},
#'  uses value determined in \code{"run_sensitivity_analysis"}.
#' @param factor_type Character; one of \code{"common_top_factor_features"}
#' or \code{"top_factor_features"}
#'   (for \code{"factor_features"}).
#' @param variable_map A data.frame with columns \code{exp_name} and
#' \code{feature} used when \code{feature_type = "omics"} to select
#'  features per experiment.
#' @param padj_method Multiple-testing correction method passed to
#' \code{stats::p.adjust} (default \code{"fdr"}).
#' @param pval_thresh Adjusted p-value filter applied to enrichment results
#' (default \code{0.1}).
#' @param min_set Minimum overlap size filter on enriched terms
#' (default \code{3}).
#' @param max_set Maximum overlap size filter on enriched terms
#' (default \code{800}).
#'
#' @return A \code{data.frame} (tibble) of enrichment results across
#' experiments with columns returned by \code{.run_fenr()} plus
#' \code{exp_name} and \code{padj}, filtered by \code{pval_thresh},
#' \code{min_set}, and \code{max_set}. If no results, returns an empty tibble.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment experiments metadata
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr filter inner_join pull mutate bind_rows select
#' @importFrom rlang sym `%||%`
#' @importFrom stats p.adjust
#' @import fenr
.run_feature_enr <- function(
    expomicset,
    feature_type,
    db,
    fdata,
    species,
    fenr_col,
    feature_col,
    deg_pval_col,
    deg_pval_threshold,
    deg_logfc_col,
    deg_logfc_threshold,
    score_col,
    score_threshold,
    factor_type,
    variable_map,
    padj_method,
    pval_thresh,
    min_set,
    max_set) {
    # Get unique experiments
    exps <- names(MultiAssayExperiment::experiments(expomicset))

    # Set variable for enrichment results
    enr_res <- list()

    for (exp in exps) {
        sprintf(
            "Working on '%s' %s %s enrichment analysis",
            feature_type, exp, db
        )

        # Subset universe from rowData of experiment
        universe_genes <- SummarizedExperiment::rowData(
            expomicset[[exp]]
        )[[feature_col]] |>
            unique()

        # Select features per method
        selected_genes <- switch(feature_type,
            "degs" = {
                expomicset@metadata |>
                    purrr::pluck("differential_analysis") |>
                    purrr::pluck("differential_abundance") |>
                    dplyr::filter(exp_name == exp) |>
                    dplyr::filter(
                        !!rlang::sym(deg_pval_col) < deg_pval_threshold,
                        abs(!!rlang::sym(deg_logfc_col)) > deg_logfc_threshold
                    ) |>
                    pull(!!rlang::sym(feature_col)) |>
                    unique()
            },
            "degs_robust" = {
                sens <- expomicset@metadata |>
                    purrr::pluck("differential_analysis") |>
                    purrr::pluck("sensitivity_analysis")

                sens |>
                    purrr::pluck("feature_stability") |>
                    dplyr::filter(
                        exp_name == exp,
                        !!rlang::sym(score_col) > (
                            score_threshold %||% sens$score_thresh)
                    ) |>
                    dplyr::select(feature, exp_name) |>
                    dplyr::inner_join(
                        fdata |>
                            dplyr::filter(.exp_name == exp),
                        by = c(
                            "feature" = ".feature",
                            "exp_name" = ".exp_name"
                        )
                    ) |>
                    dplyr::pull(!!dplyr::sym(feature_col)) |>
                    unique()
            },
            "factor_features" = {
                expomicset@metadata |>
                    purrr::pluck("differential_analysis") |>
                    purrr::pluck("multiomics_integration") |>
                    purrr::pluck(factor_type) |>
                    dplyr::filter(exp_name == exp) |>
                    dplyr::inner_join(
                        fdata |>
                            dplyr::filter(.exp_name == exp),
                        by = c(
                            "feature" = ".feature",
                            "exp_name" = ".exp_name"
                        )
                    ) |>
                    dplyr::pull(!!dplyr::sym(feature_col)) |>
                    unique()
            },
            "omics" = {
                variable_map |>
                    dplyr::filter(exp_name == exp) |>
                    dplyr::inner_join(
                        fdata |>
                            dplyr::filter(.exp_name == exp),
                        by = c(
                            "feature" = ".feature",
                            "exp_name" = ".exp_name"
                        )
                    ) |>
                    dplyr::pull(!!dplyr::sym(feature_col)) |>
                    unique()
            }
        )

        if (length(selected_genes) < 1) {
            enr <- NULL
        } else if (is.null(universe_genes) | all(is.na(universe_genes))) {
            enr <- NULL
        } else {
            # Run enrichment
            enr <- .run_fenr(
                selected_genes = selected_genes,
                universe_genes = universe_genes,
                db = db,
                species = species,
                fenr_col = fenr_col
            ) |>
                dplyr::mutate(exp_name = exp)
        }

        enr_res[[exp]] <- enr
    }

    if (all(purrr::map_lgl(enr_res, is.null))) {
        return(tibble::tibble())
    }

    # Multiple-hypothesis correction after all tests are done
    enr_res <- enr_res |>
        dplyr::bind_rows() |>
        dplyr::mutate(padj = p.adjust(p_value,
            method = padj_method
        )) |>
        dplyr::filter(padj < pval_thresh) |>
        dplyr::filter(
            n_with_sel > min_set,
            n_with_sel < max_set
        )

    return(enr_res)
}
# --- Run Enrichment on Exposure-Correlated Features ---------
#' Run correlation-based enrichment analysis for each (experiment, category)
#'
#' This internal helper performs enrichment analysis using correlation-derived
#' features across experiments and exposure categories.
#' It handles feature types ending in \code{"_cor"}:
#' \itemize{
#'   \item \code{"degs_cor"}: correlation between DEGs and exposures
#'   \item \code{"omics_cor"}: correlation between omics features and exposures
#'   \item \code{"factor_features_cor"}: correlation between factor features and exposures
#' }
#'
#' For each experiment by category pair, it:
#' \enumerate{
#'   \item Extracts the gene universe from \code{rowData()}.
#'   \item Selects features based on correlation metadata and
#'   \code{pivot_feature()}.
#'   \item Runs \code{.run_fenr()} with the selected genes.
#'   \item Annotates results with \code{exp_name} and \code{category}.
#' }
#' Then, it applies multiple testing correction and filters enriched
#' terms by adjusted p-value and gene set size.
#'
#' @inheritParams .run_feature_enr
#' @param feature_type One of \code{"degs_cor"}, \code{"omics_cor"},
#'  \code{"factor_features_cor"}.
#'
#' @return A data frame of enrichment results across experiments and
#' categories,
#' with columns returned by \code{.run_fenr()} plus \code{exp_name},
#'  \code{category}, and \code{padj}.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom purrr pluck
#' @importFrom MultiAssayExperiment metadata experiments
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr filter inner_join pull mutate bind_rows select
#' @importFrom rlang sym
#' @importFrom stats p.adjust
#' @import fenr
.run_feature_cor_enr <- function(
    expomicset,
    feature_type,
    db,
    fdata,
    species,
    fenr_col,
    feature_col,
    deg_pval_col,
    deg_pval_threshold,
    deg_logfc_col,
    deg_logfc_threshold,
    score_col,
    score_threshold,
    factor_type,
    variable_map,
    padj_method,
    pval_thresh,
    min_set,
    max_set) {
    # Set the categories and experiment names to loop through
    cor_table <- expomicset@metadata |>
        purrr::pluck("correlation") |>
        purrr::pluck(gsub("_.*", "", feature_type))

    exps <- unique(cor_table$exp_name)
    categories <- unique(cor_table$category)

    # Set variable for enrichment results
    enr_res <- list()

    for (cat_i in categories) {
        for (exp in exps) {
            sprintf(
                "Working on '%s' %s %s %s enrichment analysis",
                feature_type, exp, cat_i, db
            )
            # Subset universe from rowData of experiment
            universe_genes <- SummarizedExperiment::rowData(
                expomicset[[exp]]
            )[[feature_col]] |>
                unique()

            # Select features per feature type
            selected_genes <- switch(feature_type,
                "degs_cor" = {
                    expomicset@metadata |>
                        purrr::pluck("correlation") |>
                        purrr::pluck("degs") |>
                        dplyr::filter(
                            exp_name == exp,
                            category == cat_i
                        ) |>
                        dplyr::inner_join(
                            fdata |>
                                dplyr::filter(.exp_name == exp),
                            by = c(
                                "feature" = ".feature",
                                "exp_name" = ".exp_name"
                            )
                        ) |>
                        dplyr::pull(!!dplyr::sym(feature_col)) |>
                        unique()
                },
                "omics_cor" = {
                    expomicset@metadata |>
                        purrr::pluck("correlation") |>
                        purrr::pluck("omics") |>
                        dplyr::filter(
                            exp_name == exp,
                            category == cat_i
                        ) |>
                        dplyr::inner_join(
                            fdata |>
                                dplyr::filter(.exp_name == exp),
                            by = c(
                                "feature" = ".feature",
                                "exp_name" = ".exp_name"
                            )
                        ) |>
                        dplyr::pull(!!dplyr::sym(feature_col)) |>
                        unique()
                },
                "factor_features_cor" = {
                    expomicset@metadata |>
                        purrr::pluck("correlation") |>
                        purrr::pluck("factor_features") |>
                        dplyr::filter(
                            exp_name == exp,
                            category == cat_i
                        ) |>
                        dplyr::inner_join(
                            fdata |>
                                dplyr::filter(.exp_name == exp),
                            by = c(
                                "feature" = ".feature",
                                "exp_name" = ".exp_name"
                            )
                        ) |>
                        dplyr::pull(!!dplyr::sym(feature_col)) |>
                        unique()
                }
            )

            if (length(selected_genes) < 1) {
                enr <- NULL
            } else if (is.null(universe_genes) | all(is.na(universe_genes))) {
                enr <- NULL
            } else {
                # Run enrichment
                enr <- .run_fenr(
                    selected_genes = selected_genes,
                    universe_genes = universe_genes,
                    db = db,
                    species = species,
                    fenr_col = fenr_col
                ) |>
                    dplyr::mutate(exp_name = exp) |>
                    dplyr::mutate(category = cat_i)
            }

            enr_res[[exp]] <- enr
        }
    }

    if (all(purrr::map_lgl(enr_res, is.null))) {
        return(tibble::tibble())
    }

    # Multiple-hypothesis correction after all tests are done
    enr_res <- enr_res |>
        dplyr::bind_rows() |>
        dplyr::mutate(padj = p.adjust(p_value,
            method = padj_method
        )) |>
        dplyr::filter(padj < pval_thresh) |>
        dplyr::filter(
            n_with_sel > min_set,
            n_with_sel < max_set
        )

    return(enr_res)
}

# --- Group Enrichment Terms Function ------------
#' Group enriched GO terms into clusters based on gene overlap
#'
#' This internal function clusters enriched GO terms based on pairwise Jaccard
#' overlap of their gene lists. It is typically used downstream of enrichment
#' analysis to simplify redundant terms into groups.
#'
#' @param enr_res A data frame containing enrichment results. Must contain:
#'   - `term_name`: term identifier (e.g. GO term)
#'   - `ids`: comma-separated list of associated genes
#' @param clustering_approach Clustering method to use. Defaults to `"hclust"`.
#'
#' @return A data frame with an added `go_group` column indicating the cluster
#'   assignment for each term.
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import tibble
#' @importFrom rlang set_names
#' @keywords internal
#' @noRd
.group_enr_res <- function(enr_res, clustering_approach = "diana") {
    if (nrow(enr_res) < 2) {
        message("GO clustering skipped - fewer than 2 enriched terms.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    message("Determining Number of GO Term Clusters...")

    go_term_list <- enr_res |>
        dplyr::mutate(gene_list = stringr::str_split(ids, ",")) |>
        tidyr::unnest(gene_list) |>
        dplyr::mutate(gene_list = stringr::str_trim(gene_list)) |>
        dplyr::distinct(term_name, gene_list) |>
        dplyr::group_by(term_name) |>
        dplyr::summarise(genes = list(unique(gene_list)), .groups = "drop")

    if (nrow(go_term_list) < 2) {
        message("GO clustering skipped - fewer than 2 unique terms.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    go_groups <- go_term_list$genes |>
        rlang::set_names(go_term_list$term_name) |>
        .get_pairwise_overlaps() |>
        dplyr::select(source, target, jaccard) |>
        tidyr::pivot_wider(names_from = "target", values_from = "jaccard") |>
        tibble::column_to_rownames("source") |>
        as.matrix() |>
        .cluster_mat(clustering_approach = clustering_approach) |>
        (\(x) data.frame(
            term_name = names(x),
            go_group = paste0("Group_", as.numeric(x))
        ))()

    if (is.null(go_groups) || nrow(go_groups) == 0) {
        message("GO clustering failed - no clusters found.")
        enr_res$go_group <- rep(NA_character_, nrow(enr_res))
        return(enr_res)
    }

    dplyr::inner_join(enr_res, go_groups, by = "term_name")
}
