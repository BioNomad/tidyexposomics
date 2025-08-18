#' Run Correlation Analysis
#'
#' Computes correlations between exposures and feature types including DEGs,
#'  omics, latent factors,
#' top factor features, or principal components (PCs). Optionally computes
#' feature–feature correlations
#' to support network analysis.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param feature_type Type of features to correlate. One of `"degs"`,
#'  `"omics"`, `"factors"`, `"factor_features"`, `"exposures"`, or `"pcs"`.
#' @param exposure_cols Optional vector of exposure column names
#' (from `colData`) to use.
#' @param variable_map Optional mapping of features to include by assay
#'  for `omics` mode.
#' @param n_pcs Number of PCs to use when `feature_type = "pcs"`.
#' @param feature_cors Logical; if `TRUE`, compute correlations between
#' features rather than with exposures.
#' @param robust Logical; restrict DEGs to those passing sensitivity threshold.
#' @param score_col Column name in sensitivity analysis with feature
#' stability score.
#' @param score_thresh Threshold for filtering robust features.
#' @param correlation_method One of `"pearson"`, `"spearman"`, or `"kendall"`.
#' @param correlation_cutoff Minimum absolute correlation to retain.
#' @param cor_pval_column Column in output to filter by p-value
#' (default: `"p.value"`).
#' @param pval_cutoff Maximum p-value or FDR threshold to retain a correlation.
#' @param deg_pval_col Column with DEG adjusted p-values.
#' @param deg_logfc_col Column with DEG log fold-changes.
#' @param deg_pval_thresh P-value cutoff for DEGs.
#' @param deg_logfc_thresh Log fold-change cutoff for DEGs.
#' @param batch_size Number of features to process per batch (default: 1500).
#' @param action Whether to `"add"` results to metadata or `"get"`
#' as a data frame.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`
#'  with results added to metadata.
#'         If `action = "get"`, returns a tidy `data.frame` of correlations.
#'
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run correlation analysis
#' mae <- mae |>
#'     run_correlation(
#'         feature_type = "exposures",
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' @export
run_correlation <- function(
    expomicset,
    feature_type = c(
        "degs",
        "omics",
        "factors",
        "factor_features",
        "exposures",
        "pcs"
    ),
    exposure_cols = NULL,
    variable_map = NULL,
    n_pcs = NULL,
    feature_cors = FALSE,
    robust = FALSE,
    score_col = "stability_score",
    score_thresh = NULL,
    correlation_method = "spearman",
    correlation_cutoff = 0.3,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05,
    deg_pval_col = "adj.P.Val",
    deg_logfc_col = "logFC",
    deg_pval_thresh = 0.05,
    deg_logfc_thresh = log2(1.5),
    batch_size = 1500,
    action = c("add", "get")) {
    feature_type <- match.arg(feature_type)
    action <- match.arg(action)

    col_df <- MultiAssayExperiment::colData(expomicset) |>
        as.data.frame()
    # dplyr::select(-dplyr::starts_with("PC"))

    exposures <- col_df |>
        dplyr::select(where(is.numeric))

    if (!is.null(exposure_cols)) {
        exposures <- exposures[, intersect(
            colnames(exposures),
            exposure_cols
        ),
        drop = FALSE
        ]
    }
    if (ncol(exposures) == 0) stop("No numeric exposures found.")

    feature_matrix <- switch(feature_type,
        degs = .extract_deg_matrix(
            expomicset,
            robust,
            score_col,
            score_thresh,
            deg_pval_col,
            deg_logfc_col,
            deg_pval_thresh,
            deg_logfc_thresh
        ),
        omics = .extract_omics_matrix(
            expomicset,
            variable_map
        ),
        factors = .extract_factor_matrix(expomicset),
        exposures = .extract_exposure_matrix(
            col_df,
            exposure_cols
        ),
        factor_features = .extract_factor_feature_matrix(expomicset),
        pcs = .extract_pc_matrix(col_df, n_pcs = n_pcs)
    )

    exposures <- exposures |>
        tibble::rownames_to_column("id")

    if (feature_type == "exposures") {
        merged_data <- exposures
        exposure_vars <- setdiff(colnames(merged_data), "id")
        feature_vars <- exposure_vars
    } else {
        if (!"id" %in% colnames(feature_matrix)) {
            feature_matrix <- feature_matrix |>
                tibble::rownames_to_column("id")
        }
        merged_data <- dplyr::left_join(exposures, feature_matrix, by = "id") |>
            na.omit()
        exposure_vars <- setdiff(colnames(exposures), "id")
        feature_vars <- setdiff(colnames(feature_matrix), "id")
    }


    correlation_df <- .run_and_clean_correlation(
        expomicset = expomicset,
        merged_data = merged_data,
        feature_vars = feature_vars,
        exposure_vars = exposure_vars,
        feature_type = feature_type,
        feature_cors = feature_cors,
        correlation_method = correlation_method,
        correlation_cutoff = correlation_cutoff,
        cor_pval_column = cor_pval_column,
        pval_cutoff = pval_cutoff,
        batch_size = batch_size
    )

    if (action == "add") {
        if (feature_cors) {
            all_metadata <- MultiAssayExperiment::metadata(expomicset)
            res_name <- paste0(feature_type, "_feature_cor")
            all_metadata$correlation[[res_name]] <- correlation_df
            MultiAssayExperiment::metadata(expomicset) <- all_metadata

            step_record <- list(
                run_correlation = list(
                    timestamp = Sys.time(),
                    params = list(
                        feature_type = feature_type,
                        correlation_method = correlation_method,
                        correlation_cutoff = correlation_cutoff,
                        pval_cutoff = pval_cutoff
                    ),
                    notes = paste0(
                        "Correlated ",
                        feature_type,
                        " features with features."
                    )
                )
            )
            names(step_record) <- paste0("run_correlation_", feature_type)
        } else {
            all_metadata <- MultiAssayExperiment::metadata(expomicset)
            all_metadata$correlation[[feature_type]] <- correlation_df
            MultiAssayExperiment::metadata(expomicset) <- all_metadata

            step_record <- list(
                run_correlation = list(
                    timestamp = Sys.time(),
                    params = list(
                        feature_type = feature_type,
                        correlation_method = correlation_method,
                        correlation_cutoff = correlation_cutoff,
                        pval_cutoff = pval_cutoff
                    ),
                    notes = paste0(
                        "Correlated ",
                        feature_type,
                        " features with exposures."
                    )
                )
            )
            names(step_record) <- paste0("run_correlation_", feature_type)
        }


        MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(expomicset)$summary$steps,
            step_record
        )
        return(expomicset)
    } else {
        return(correlation_df)
    }
}

# --- Run exposure–feature correlation function ----------
#' Compute exposure–feature correlations in batches
#'
#' Internal helper that correlates a set of exposure variables against a set
#' of feature variables by processing features in batches. Uses
#' \code{Hmisc::rcorr()} to compute correlation coefficients and p-values, then
#' filters by absolute correlation magnitude and a p-value/FDR column.
#'
#' @param merged_data A data frame/matrix with one row per sample and columns
#'   for all \code{exposure_vars} and \code{feature_vars}. Must be numeric for
#'   selected columns. Rows with non-finite values across the selected columns
#'   are removed prior to correlation.
#' @param exposure_vars Character vector of exposure column names (subset of
#'   \code{colnames(merged_data)}).
#' @param feature_vars Character vector of feature column names (subset of
#'   \code{colnames(merged_data)}).
#' @param correlation_method Correlation method passed to \code{Hmisc::rcorr()},
#'   typically \code{"pearson"} or \code{"spearman"}.
#' @param correlation_cutoff Numeric; retain pairs with
#'   \eqn{|r| >} \code{correlation_cutoff}.
#' @param cor_pval_column Name of the p-value column to filter on after FDR
#'   computation; usually \code{"p.value"} or \code{"FDR"}.
#' @param pval_cutoff Numeric; maximum allowed value for \code{cor_pval_column}.
#' @param batch_size Integer; number of feature columns processed per batch.
#'
#' @return A \code{data.frame} with one row per retained exposure–feature pair
#' and columns: \code{var1} (exposure), \code{var2} (feature),
#' \code{correlation}, \code{p.value}, and \code{FDR}.
#' Returns an empty data frame when no pairs pass filters.
#'
#' @details
#' Correlations are computed within each batch, then results are row-bound.
#' The FDR is computed within each batch; if you require global FDR across all
#' pairs, compute \code{p.adjust} on the combined result after binding rows.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Hmisc rcorr
#' @importFrom dplyr inner_join filter mutate bind_rows
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom stats p.adjust
.run_correlation_batches <- function(
    merged_data,
    exposure_vars,
    feature_vars,
    correlation_method,
    correlation_cutoff,
    cor_pval_column,
    pval_cutoff,
    batch_size) {
    correlation_results <- list()
    batches <- split(feature_vars, ceiling(seq_along(feature_vars) / batch_size))

    for (i in seq_along(batches)) {
        features <- batches[[i]]
        cols_to_use <- intersect(c(exposure_vars, features), colnames(merged_data))
        mat <- as.matrix(merged_data[, cols_to_use, drop = FALSE])
        mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
        if (nrow(mat) < 3) next

        corr_mat <- Hmisc::rcorr(mat, type = correlation_method)

        corr_df <- as.data.frame(corr_mat$r) |>
            tibble::rownames_to_column("var1") |>
            tidyr::pivot_longer(-var1,
                names_to = "var2",
                values_to = "correlation"
            )
        pval_df <- as.data.frame(corr_mat$P) |>
            tibble::rownames_to_column("var1") |>
            tidyr::pivot_longer(-var1,
                names_to = "var2",
                values_to = "p.value"
            )

        merged <- dplyr::inner_join(corr_df,
            pval_df,
            by = c("var1", "var2")
        ) |>
            dplyr::filter(
                var1 %in% exposure_vars,
                var2 %in% features,
                abs(correlation) > correlation_cutoff
            )

        correlation_results[[i]] <- merged
    }

    dplyr::bind_rows(correlation_results) |>
        dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
        dplyr::filter(!!rlang::sym(cor_pval_column) < pval_cutoff)
}


# --- Run feature–feature correlation function ------------
#' Compute feature–feature correlations in batches
#'
#' Internal helper that correlates a (potentially large) set of feature
#' variables against itself by processing features in batches. Uses
#' \code{Hmisc::rcorr()} to compute correlation coefficients and p-values,
#' then filters by absolute correlation magnitude and an adjusted p-value
#' (FDR) or raw p-value column.
#'
#' @param feature_data A data frame/matrix with one row per sample and columns
#'   for all \code{feature_vars}. Selected columns must be numeric. Rows with
#'   non-finite values across selected columns are removed prior to
#'   correlation.
#' @param feature_vars Character vector of feature column names (subset of
#'   \code{colnames(feature_data)}).
#' @param correlation_method Correlation method passed to \code{Hmisc::rcorr()},
#'   typically \code{"pearson"} or \code{"spearman"}.
#' @param correlation_cutoff Numeric; retain pairs with
#'   \eqn{|r| >} \code{correlation_cutoff}.
#' @param cor_pval_column Name of the p-value column to filter on after FDR
#'   computation; usually \code{"p.value"} or \code{"FDR"}.
#' @param pval_cutoff Numeric; maximum allowed value for \code{cor_pval_column}.
#' @param batch_size Integer; number of feature columns processed per batch.
#'
#' @return A \code{data.frame} with one row per retained feature–feature pair
#'   and columns: \code{var1}, \code{var2}, \code{correlation}, \code{p.value},
#'   and \code{FDR}. Returns an empty data frame when no pairs pass filters.
#'
#' @details
#' Correlations are computed within each batch, then results are row-bound.
#' The FDR is computed within each batch; if you require global FDR across all
#' pairs, recompute \code{p.adjust} on the combined result after binding.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Hmisc rcorr
#' @importFrom dplyr inner_join filter mutate bind_rows
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom stats p.adjust
.run_correlation_omics_batches <- function(
    feature_data,
    feature_vars,
    correlation_method,
    correlation_cutoff,
    cor_pval_column,
    pval_cutoff,
    batch_size) {
    correlation_results <- list()
    batches <- split(feature_vars, ceiling(seq_along(feature_vars) / batch_size))

    for (i in seq_along(batches)) {
        features <- batches[[i]]
        cols_to_use <- intersect(features, colnames(feature_data))
        mat <- as.matrix(feature_data[, cols_to_use, drop = FALSE])
        mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
        if (nrow(mat) < 3) next

        corr_mat <- Hmisc::rcorr(mat, type = correlation_method)

        corr_df <- as.data.frame(corr_mat$r) |>
            tibble::rownames_to_column("var1") |>
            tidyr::pivot_longer(-var1,
                names_to = "var2",
                values_to = "correlation"
            )
        pval_df <- as.data.frame(corr_mat$P) |>
            tibble::rownames_to_column("var1") |>
            tidyr::pivot_longer(-var1,
                names_to = "var2",
                values_to = "p.value"
            )

        merged <- dplyr::inner_join(corr_df,
            pval_df,
            by = c("var1", "var2")
        ) |>
            dplyr::filter(abs(correlation) > correlation_cutoff)

        correlation_results[[i]] <- merged
    }

    dplyr::bind_rows(correlation_results) |>
        dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
        dplyr::filter(!!rlang::sym(cor_pval_column) < pval_cutoff)
}

# --- Run and Clean Correlation Function --------------
#' Run correlation and normalize variable names
#'
#' Internal helper that dispatches to the appropriate correlation routine
#' (exposure–feature vs feature–feature), then normalizes variable names by
#' stripping assay/experiment prefixes from feature columns. When correlating
#' exposures with features, the function also joins the codebook metadata.
#'
#' The output keeps \code{var1}/\code{var2} for downstream compatibility and,
#' when \code{feature_cors = FALSE}, adds convenience columns
#' \code{exposure} (== \code{var1}) and \code{feature} (== \code{var2}).
#' It also annotates \code{var1_type} and \code{var2_type}.
#'
#' @param expomicset A \code{MultiAssayExperiment}.
#' @param merged_data Data frame used for correlation (rows = samples).
#' @param feature_vars Character vector of feature column names.
#' @param exposure_vars Character vector of exposure column names.
#' @param feature_type One of \code{"degs"}, \code{"omics"}, \code{"factors"},
#'   \code{"factor_features"}, \code{"exposures"}, or \code{"pcs"}.
#' @param feature_cors Logical; if \code{TRUE}, do feature–feature correlation.
#' @param correlation_method Correlation method for \code{Hmisc::rcorr()}.
#' @param correlation_cutoff Keep pairs with \eqn{|r| >} this threshold.
#' @param cor_pval_column Name of p-value column to filter (e.g., "p.value",
#'   "FDR").
#' @param pval_cutoff Numeric cutoff applied to \code{cor_pval_column}.
#' @param batch_size Integer batch size for correlation.
#'
#' @return A data frame with \code{var1}, \code{var2}, \code{correlation},
#'   \code{p.value}, \code{FDR}, and, when applicable, \code{exposure},
#'   \code{feature}, \code{exp_name}. Includes \code{var1_type} and
#'   \code{var2_type}.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment experiments metadata
#' @importFrom dplyr left_join mutate
#' @importFrom stringr str_extract str_remove str_replace_all
.run_and_clean_correlation <- function(
    expomicset,
    merged_data,
    feature_vars,
    exposure_vars,
    feature_type,
    feature_cors,
    correlation_method,
    correlation_cutoff,
    cor_pval_column,
    pval_cutoff,
    batch_size) {
    if (feature_cors && feature_type %in% c(
        "degs",
        "omics",
        "factors",
        "factor_features"
    )) {
        # Run correlation on features
        correlation_df <- .run_correlation_omics_batches(
            merged_data,
            feature_vars,
            correlation_method,
            correlation_cutoff,
            cor_pval_column,
            pval_cutoff,
            batch_size
        )

        # Grab experiment names
        exp_names <- names(MultiAssayExperiment::experiments(expomicset))

        # Change names of the columns to match column type
        correlation_df <- correlation_df |>
            # Add in exp_name for the assay
            dplyr::mutate(
                exp_name_1 = stringr::str_extract(
                    var1,
                    paste0("^(", paste0(
                        stringr::str_replace_all(
                            exp_names, " ", "[ _]"
                        ),
                        collapse = "|"
                    ), ")")
                ),
                # Remove matched exp_name and any underscore or space after it
                var1 = dplyr::case_when(
                    !is.na(exp_name_1) ~ stringr::str_remove(
                        var1, paste0("^", exp_name_1, "[ _]?")
                    ),
                    TRUE ~ var1
                ),
                exp_name_2 = stringr::str_extract(
                    var2,
                    paste0(
                        "^(",
                        paste0(
                            stringr::str_replace_all(
                                exp_names, " ", "[ _]"
                            ),
                            collapse = "|"
                        ), ")"
                    )
                ),
                # Remove matched exp_name and any underscore or space after it
                var2 = dplyr::case_when(
                    !is.na(exp_name_2) ~ stringr::str_remove(
                        var2, paste0("^", exp_name_2, "[ _]?")
                    ),
                    TRUE ~ var2
                )
            )
    } else {
        correlation_df <- .run_correlation_batches(
            merged_data,
            exposure_vars,
            feature_vars,
            correlation_method,
            correlation_cutoff,
            cor_pval_column,
            pval_cutoff,
            batch_size
        )
    }


    if (!feature_cors && feature_type %in% c(
        "degs",
        "omics",
        "factors",
        "factor_features"
    )) {
        # Grab experiment names
        exp_names <- names(MultiAssayExperiment::experiments(expomicset))

        # Change names of the columns to match column type
        correlation_df <- correlation_df |>
            dplyr::rename(
                exposure = var1,
                feature = var2
            ) |>
            # Add in exp_name for the assay
            dplyr::mutate(
                exp_name = stringr::str_extract(
                    feature,
                    paste0(
                        "^(",
                        paste0(stringr::str_replace_all(exp_names, " ", "[ _]"),
                            collapse = "|"
                        ), ")"
                    )
                ),
                # Remove matched exp_name and any underscore or space after it
                feature = dplyr::case_when(
                    !is.na(exp_name) ~ stringr::str_remove(
                        feature, paste0("^", exp_name, "[ _]?")
                    ),
                    TRUE ~ feature
                )
            )


        # Merge with exposure metadata
        correlation_df <- correlation_df |>
            dplyr::left_join(
                MultiAssayExperiment::metadata(expomicset)$codebook,
                by = c("exposure" = "variable")
            )
    }
    return(correlation_df)
}

# --- Extract DEGs to Correlate Function -----------
#' Extract DEG expression matrix from a MultiAssayExperiment
#'
#' Internal helper that extracts log2-normalized expression values for
#' differentially expressed genes (DEGs) across assays in a
#' \code{MultiAssayExperiment}, optionally restricted to robust features
#' based on sensitivity analysis.
#'
#' @param expomicset A \code{MultiAssayExperiment} containing DEG and (optional)
#'   sensitivity metadata.
#' @param robust Logical; if \code{TRUE}, filter to robust features using
#'   sensitivity scores.
#' @param score_col Name of the column with stability scores.
#' @param score_thresh Threshold for stability score filtering.
#' @param deg_pval_col Name of column containing DEG p-values.
#' @param deg_logfc_col Name of column containing DEG log fold-changes.
#' @param deg_pval_thresh P-value threshold for DEG filtering.
#' @param deg_logfc_thresh Absolute logFC threshold for DEG filtering.
#'
#' @return A data frame of log2-normalized expression values for DEGs
#'   (columns) with one row per sample and a column \code{"id"} for sample ID.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom purrr pluck
#' @importFrom dplyr filter pull semi_join
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
.extract_deg_matrix <- function(
    expomicset,
    robust,
    score_col,
    score_thresh,
    deg_pval_col,
    deg_logfc_col,
    deg_pval_thresh,
    deg_logfc_thresh) {
    da <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "differential_analysis",
            "differential_abundance"
        )
    if (is.null(da)) stop("No differential_abundance in metadata")

    da <- da |>
        dplyr::filter(
            !!rlang::sym(deg_pval_col) < deg_pval_thresh,
            abs(!!rlang::sym(deg_logfc_col)) > deg_logfc_thresh
        )

    if (robust) {
        sens <- MultiAssayExperiment::metadata(expomicset) |>
            purrr::pluck(
                "differential_analysis",
                "sensitivity_analysis"
            )
        if (is.null(sens)) stop("No sensitivity_analysis in metadata")

        stable_feats <- sens$feature_stability |>
            dplyr::filter(
                !!dplyr::sym(score_col) > (score_thresh %||% sens$score_thresh)
            )
        da <- dplyr::semi_join(da,
            stable_feats,
            by = c("exp_name", "feature")
        )
    }

    mats <- lapply(unique(da$exp_name), function(name) {
        se <- .update_assay_colData(expomicset, name)
        feats <- da |>
            dplyr::filter(exp_name == name) |>
            dplyr::pull(feature)
        se <- se[rownames(se) %in% feats, , drop = FALSE]
        mat <- SummarizedExperiment::assay(se) |>
            t() |>
            as.data.frame()
        colnames(mat) <- paste(name, colnames(mat), sep = "_")
        tibble::rownames_to_column(mat, "id")
    })

    purrr::reduce(mats, dplyr::full_join, by = "id")
}


# --- Extract User-specified Omics to Correlate Function ---------
#' Extract omics expression matrix from MultiAssayExperiment
#'
#' Internal helper that extracts log2-transformed expression matrices for
#' omics features in a \code{MultiAssayExperiment}, optionally restricted to
#' a subset defined by a variable mapping table.
#'
#' @param expomicset A \code{MultiAssayExperiment} object with omics assays.
#' @param variable_map Optional data frame with columns \code{variable} and
#'   \code{exp_name} indicating which features to extract from which assays.
#'
#' @return A data frame with one row per sample and one column per omics
#'   feature, prefixed by assay name. Includes a column \code{"id"} for sample ID.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment experiments
#' @importFrom purrr imap compact reduce
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
.extract_omics_matrix <- function(expomicset, variable_map) {
    log2_omics <- .log2_multiassay(expomicset)
    selected <- if (!is.null(variable_map)) {
        split(
            variable_map$variable,
            variable_map$exp_name
        )
    } else {
        lapply(MultiAssayExperiment::experiments(log2_omics), rownames)
    }

    dfs <- purrr::imap(
        MultiAssayExperiment::experiments(log2_omics),
        function(se, name) {
            feats <- selected[[name]]
            if (is.null(feats)) {
                return(NULL)
            }
            se <- se[rownames(se) %in% feats, , drop = FALSE]
            df <- SummarizedExperiment::assay(se) |>
                t() |>
                as.data.frame()
            colnames(df) <- paste(name, colnames(df), sep = "_")
            tibble::rownames_to_column(df, "id")
        }
    ) |> purrr::compact()

    purrr::reduce(dfs,
        dplyr::full_join,
        by = "id"
    )
}

# --- Extract Latent Factor to Correlate Function ---------
#' Extract sample × factor matrix from multi-omics integration
#'
#' Internal helper that retrieves sample-level factor scores from the stored
#' multi-omics integration results in a \code{MultiAssayExperiment}. Supports
#' \code{MOFA}, \code{MCIA}, and \code{MCCA}. Returns a data frame with one row
#' per sample and a column \code{"id"} for sample IDs.
#'
#' @param expomicset A \code{MultiAssayExperiment} with integration results
#'   saved in \code{metadata(expomicset)$multiomics_integration$integration_results}.
#'
#' @return A data frame of sample × factor scores with an \code{"id"} column.
#'
#' @details
#' For MOFA, factors are obtained via \code{MOFA2::get_factors()} (first view).
#' For MCIA, the \code{@global_scores} slot is used. For MCCA, the
#' \code{@sample_scores} slot is used. An error is thrown for unsupported
#' methods.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom purrr pluck
#' @importFrom tibble rownames_to_column
.extract_factor_matrix <- function(expomicset) {
    result <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "multiomics_integration",
            "integration_results"
        )
    mat <- if (result$method == "MOFA") {
        .check_suggested("MOFA2")
        MOFA2::get_factors(result$result)[[1]]
    } else if (result$method == "MCIA") {
        result$result@global_scores
    } else if (result$method == "MCCA") {
        result$result@sample_scores
    } else {
        stop("Unsupported integration method")
    }
    as.data.frame(mat) |> tibble::rownames_to_column("id")
}

# --- Extract Latent Factor Features to Correlate Function ---------
#' Extract matrix for top factor features across assays
#'
#' Internal helper that builds a sample × feature matrix for features selected
#' from multi-omics integration (stored in
#' \code{metadata(expomicset)$multiomics_integration$top_factor_features}).
#' For each assay, it subsets to the specified features, extracts log2 values,
#' and returns a joined data frame with an \code{"id"} column.
#'
#' @param expomicset A \code{MultiAssayExperiment} containing omics assays and
#'   multi-omics integration metadata with a \code{top_factor_features} table
#'   that includes \code{exp_name} and \code{feature} columns.
#'
#' @return A data frame with one row per sample and columns for selected
#'   features (prefixed by assay name), plus an \code{"id"} column.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom MultiAssayExperiment metadata experiments
#' @importFrom purrr pluck imap compact reduce
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
.extract_factor_feature_matrix <- function(expomicset) {
    top_feats <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "multiomics_integration",
            "top_factor_features"
        )

    if (is.null(top_feats)) stop("No top_factor_features found in metadata.")
    if (!all(c("exp_name", "feature") %in% colnames(top_feats))) {
        stop("top_factor_features must contain 'exp_name' and 'feature' columns.")
    }

    log2_omics <- .log2_multiassay(expomicset)
    selected <- split(top_feats$feature, top_feats$exp_name)

    dfs <- purrr::imap(
        MultiAssayExperiment::experiments(log2_omics),
        function(se, name) {
            feats <- selected[[name]]
            if (is.null(feats)) {
                return(NULL)
            }
            se <- se[rownames(se) %in% feats, , drop = FALSE]
            df <- SummarizedExperiment::assay(se) |>
                t() |>
                as.data.frame()
            colnames(df) <- paste(name, colnames(df), sep = "_")
            tibble::rownames_to_column(df, "id")
        }
    ) |> purrr::compact()

    purrr::reduce(dfs, dplyr::full_join, by = "id")
}

# --- Extract Exposure Data to Correlate Function -----------
#' Extract numeric exposures from colData
#'
#' Internal helper that selects numeric exposure variables from a sample-level
#' data frame (e.g., \code{colData(expomicset)}). Optionally restricts to a
#' user-provided set of exposure column names.
#'
#' @param col_df A data frame of sample metadata (rows = samples).
#' @param exposure_cols Optional character vector of exposure column names to
#'   retain. If \code{NULL}, all numeric columns are used.
#'
#' @return A data frame with one row per sample containing the selected numeric
#'   exposure columns and an \code{"id"} column of sample IDs.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr select where
#' @importFrom tibble rownames_to_column
.extract_exposure_matrix <- function(col_df, exposure_cols) {
    df <- col_df |> dplyr::select(where(is.numeric))
    if (!is.null(exposure_cols)) df <- df[, exposure_cols, drop = FALSE]
    tibble::rownames_to_column(df, "id")
}

# --- Extract Principal Components to Correlate Function ---------
#' Extract principal component (PC) matrix from metadata
#'
#' Internal helper that selects numeric columns named like \code{PC1}, \code{PC2},
#' etc., from a sample-level data frame (e.g., \code{colData}). Optionally
#' restricts to the first \code{n_pcs} PCs by numeric index.
#'
#' @param col_df A data frame of sample metadata (rows = samples).
#' @param n_pcs Optional integer. If provided, keep only the first \code{n_pcs}
#'   PCs by their numeric suffix.
#'
#' @return A data frame with one row per sample containing selected PC columns
#'   and an \code{"id"} column of sample IDs.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr select where matches
#' @importFrom readr parse_number
#' @importFrom tibble rownames_to_column
.extract_pc_matrix <- function(col_df, n_pcs = NULL) {
    all_pcs <- col_df |>
        dplyr::select(where(is.numeric)) |>
        dplyr::select(matches("^PC\\d+$"))

    if (!is.null(n_pcs)) {
        pc_names <- colnames(all_pcs)
        # Sort by numeric PC number
        pc_names_sorted <- pc_names[order(readr::parse_number(pc_names))]
        selected <- pc_names_sorted[seq_len(min(n_pcs, length(pc_names_sorted)))]
        all_pcs <- all_pcs[, selected, drop = FALSE]
    }

    tibble::rownames_to_column(all_pcs, "id")
}
