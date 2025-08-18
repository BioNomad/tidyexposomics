#' Run Association Analysis
#'
#' Perform GLM-based association testing between a specified outcome and
#' features from exposures, omics,
#' latent factors, or GO PCs. Automatically adjusts for covariates and
#' supports both Gaussian and binomial models.
#'
#' @param expomicset A `MultiAssayExperiment` object containing data
#' and metadata.
#' @param outcome The outcome variable name (must be in `colData`).
#' @param source Source of features to test. One of `"omics"`,
#' `"exposures"`, `"factors"`.
#' @param covariates Optional vector of covariate names to include in the model.
#' @param feature_set Optional character vector of exposure or GO terms to test.
#' @param log_trans Optional boolean value dictating whether or not to log2
#' transform omics features.
#' @param top_n Optional integer: if using omics source, select top `n`
#'  most variable features.
#' @param family GLM family; `"gaussian"` or `"binomial"`.
#' @param correction_method Method for p-value adjustment (default: `"fdr"`).
#' @param action If `"add"` (default), saves results to metadata;
#' else returns results as list.
#' @param min_genes Minimum number of genes required to compute GO PCs.
#' @param feature_col If using GO PCs, the column in `rowData`
#' for matching gene symbols or IDs.
#' @param mirna_assays Optional character vector of assays to exclude
#' when extracting GO terms.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`.
#'  Otherwise, returns a list of:
#' - `results_df`: tidy summary of associations
#' - `covariates`: the covariates used
#' - `model_data`: model matrix used in the GLMs
#'
#'
#' @examples
#'
#' #' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run association tests
#' mae <- mae |>
#'     run_association(
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         outcome = "smoker",
#'         covariates = c("age"),
#'         family = "binomial"
#'     )
#' @importFrom stats binomial glm logLik as.formula coef family
#' @importFrom rlang .data
#' @export
run_association <- function(
    expomicset,
    outcome,
    source = c("omics", "exposures", "factors"),
    covariates = NULL,
    feature_set = NULL,
    log_trans = TRUE,
    top_n = NULL,
    family = "gaussian",
    correction_method = "fdr",
    action = "add",
    min_genes = 10,
    feature_col = NULL,
    mirna_assays = NULL) {
    # Validate inputs
    source <- match.arg(source)

    # grab coldata and scale if numeric
    data <- expomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame() |>
        dplyr::mutate_if(is.numeric, ~ as.numeric(scale(.)))

    # switch based on input
    features_df <- switch(source,
        omics = .extract_omics_features(
            expomicset,
            log_trans = log_trans,
            top_n
        ),
        exposures = .extract_exposures(
            data,
            feature_set
        ),
        factors = .extract_latent_factors(expomicset)
        # go_pcs = .extract_go_pcs(expomicset,
        #     geneset = feature_set,
        #     covariates,
        #     min_genes = min_genes,
        #     feature_col = feature_col,
        #     mirna_assays = mirna_assays
        # )
    )

    # create the model data
    if (source == "exposures") {
        model_data <- data
        feature_cols <- feature_set
    } else {
        model_data <- dplyr::left_join(
            tibble::rownames_to_column(data, "id"),
            features_df,
            by = "id"
        ) |>
            tibble::column_to_rownames("id")

        feature_cols <- setdiff(colnames(features_df), "id")
    }


    # ensure outcome variable is factor if binomial
    if (family == "binomial") {
        model_data[[outcome]] <- as.factor(model_data[[outcome]])

        if (length(unique(model_data[[outcome]])) != 2) {
            stop("Binary outcome must have exactly 2 levels")
        }
    }

    # ensure that covariates are have more than one value
    if (!is.null(covariates)) {
        covariates <- covariates[covariates %in% colnames(model_data)]
        if (length(covariates) == 0) {
            stop("No valid covariates provided.")
        }
        for (cov in covariates) {
            if (length(unique(model_data[[cov]])) <= 1) {
                # stop(paste("Covariate", cov, "must have more than one unique value."))
                stop(sprintf("Covariate %s must have more than one unique value.", cov))
            }
        }
    }

    message("Running GLMs.")

    feature_cols <- setdiff(colnames(features_df), "id")

    results <- purrr::map_dfr(feature_cols, function(fcol) {
        fmla <- as.formula(
            paste(
                outcome, "~", fcol,
                if (!is.null(covariates)) {
                    paste("+", paste(covariates, collapse = "+"))
                } else {
                    ""
                }
            )
        )

        model <- tryCatch(
            glm(fmla, data = model_data, family = family),
            error = function(e) NULL
        )

        # Only continue if model is valid
        if (is.null(model)) {
            return(NULL)
        }

        # extract results
        model_summary <- broom::tidy(model)
        model_filtered <- model_summary |>
            dplyr::filter(term == fcol)

        if (nrow(model_filtered) == 0) {
            message(sprintf(
                "Skipping %s - term not found in model (possibly dropped)",
                fcol
            ))
            return(NULL)
        }


        if (source %in% c("exposures", "factors")) {

        }
        # Compute R^2
        model_filtered$r2 <- .calc_r2(
            model,
            model_data,
            family = family,
            outcome = outcome
        )$r2
        model_filtered$adj_r2 <- .calc_r2(
            model,
            model_data,
            family = family,
            outcome = outcome
        )$adj_r2
        return(model_filtered)
    }) |>
        dplyr::mutate(
            p_adjust = p.adjust(
                p.value,
                method = correction_method
            )
        ) |>
        dplyr::mutate(outcome = outcome)

    # Annotate the results based on the source
    results <- .annotate_results_by_source(
        results = results,
        source = source,
        expomicset = expomicset
    )

    if (
        (paste0("assoc_", source) %in%
            names(MultiAssayExperiment::metadata(expomicset)$association)) &&
            identical(
                covariates,
                MultiAssayExperiment::metadata(expomicset)$association[[
                    paste0("assoc_", source)
                ]]$covariates
            )) {
        if (any(results$term %in%
            MultiAssayExperiment::metadata(expomicset)$association[[
                paste0("assoc_", source)
            ]]$results_df$term)) {
            stop("Association results for this feature already exists.")
        } else {
            results <- expomicset |>
                MultiAssayExperiment::metadata() |>
                purrr::pluck("association") |>
                purrr::pluck(paste0("assoc_", source)) |>
                purrr::pluck("results_df") |>
                bind_rows(results)
        }
    }

    # add results to metadata
    if (action == "add") {
        MultiAssayExperiment::metadata(expomicset)$association[[
            paste0("assoc_", source)
        ]] <- list(
            results_df = results,
            covariates = covariates
        )

        # Add in a step record
        step_record <- list(
            run_association = list(
                timestamp = Sys.time(),
                params = list(
                    outcome = outcome,
                    source = source,
                    covariates = covariates,
                    feature_set = feature_set,
                    top_n = top_n,
                    family = family,
                    correction_method = correction_method,
                    min_genes = min_genes,
                    feature_col = feature_col,
                    mirna_assays = mirna_assays
                ),
                notes = paste0(
                    "Performed association analysis using source: ",
                    source
                )
            )
        )

        MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(expomicset)$summary$steps,
            step_record
        )

        return(expomicset)
    } else {
        return(list(
            results_df = results,
            covariates = covariates,
            model_data = model_data
        ))
    }
}

# --- Extract Omics Features Function ----------
#' Extract Omics Features
#'
#' Helper for extracting and optionally selecting top omics features from a
#'  MultiAssayExperiment
#' for modeling.
#'
#' @keywords internal
#' @noRd
.extract_omics_features <- function(expomicset,
                                    log_trans = TRUE,
                                    top_n = NULL) {
    if (log_trans) {
        log2_assays <- .log2_multiassay(expomicset)
    } else {
        log2_assays <- expomicset
    }

    selected <- if (!is.null(top_n)) {
        .top_var_multiassay(log2_assays, n = top_n)
    } else {
        lapply(MultiAssayExperiment::experiments(log2_assays), rownames)
    }

    scaled <- .scale_multiassay(log2_assays,
        log2 = FALSE
    )

    dfs <- purrr::imap(
        MultiAssayExperiment::experiments(scaled),
        function(se, name) {
            df <- SummarizedExperiment::assay(
                se[selected[[name]], , drop = FALSE]
            ) |>
                t() |>
                as.data.frame()

            names(df) <- paste0(name, "_", names(df)) |>
                (\(chr) gsub(" |-", "_", chr))()

            tibble::rownames_to_column(df, "id")
        }
    )

    purrr::reduce(dfs, dplyr::full_join, by = "id")
}

# --- Extract Exposures Function ----------
#' Extract Exposure Features
#'
#' Subset exposure features (from colData) for modeling.
#'
#' @keywords internal
#' @noRd
.extract_exposures <- function(data, feature_set) {
    df <- data[, feature_set, drop = FALSE] |> as.data.frame()
    tibble::rownames_to_column(df, "id")
}

# --- Extract Latent Factors Function -------------
#' Extract Latent Factors
#'
#' Retrieve latent factors (MOFA, MCIA, MCCA, DIABLO, or RGCCA) for
#' modeling from metadata.
#'
#' @keywords internal
#' @noRd
.extract_latent_factors <- function(expomicset) {
    result <- expomicset |>
        MultiAssayExperiment::metadata() |>
        purrr::pluck("multiomics_integration") |>
        purrr::pluck("integration_results")

    mat <- switch(result$method,
        "MOFA" = {
            .check_mofa_safe()
            MOFA2::get_factors(result$result)[[1]]
        },
        "MCIA" = result$result@global_scores,
        "MCCA" = {
            result$result$sample_scores |>
                purrr::map(~ .x |>
                    as.data.frame() |>
                    tibble::rownames_to_column("sample") |>
                    tidyr::pivot_longer(-sample,
                        names_to = "factor",
                        values_to = "weight"
                    )) |>
                dplyr::bind_rows(.id = "exp_name") |>
                dplyr::group_by(factor, sample) |>
                dplyr::reframe(weight = mean(weight), .groups = "drop") |>
                tidyr::pivot_wider(names_from = "factor", values_from = "weight") |>
                tibble::column_to_rownames("sample") |>
                (\(df) {
                    colnames(df) <- gsub(" ", "_", colnames(df))
                    df
                })()
        },
        "DIABLO" = {
            result$result$variates |>
                (\(lst) lst[names(lst) != "Y"])() |>
                purrr::imap_dfc(~ {
                    comp_names <- paste(.y, colnames(.x), sep = " ")
                    df <- as.data.frame(.x)
                    colnames(df) <- comp_names
                    df
                }) |>
                tibble::rownames_to_column("id") |>
                tibble::column_to_rownames("id") |>
                (\(df) {
                    colnames(df) <- gsub(" ", "_", colnames(df))
                    df
                })()
        },
        "RGCCA" = {
            purrr::imap_dfc(result$result$Y, ~ {
                comp_names <- paste(.y, colnames(.x), sep = " ")
                df <- as.data.frame(.x)
                colnames(df) <- comp_names
                df
            }) |>
                tibble::rownames_to_column("id") |>
                tibble::column_to_rownames("id") |>
                (\(df) {
                    colnames(df) <- gsub(" ", "_", colnames(df))
                    df
                })()
        },
        stop("Unsupported integration method.")
    )

    mat <- as.data.frame(mat)
    tibble::rownames_to_column(mat, "id")
}


# --- Extract GO Principal Components Function -------
#' Extract GO Principal Components
#'
#' For each GO group in the enrichment results, compute
#' the first PC of the associated genes.
#'
# .extract_go_pcs <- function(
#     expomicset, geneset,
#     covariates,
#     min_genes = 10,
#     feature_col = NULL,
#     mirna_assays = NULL) {
#     enrich_res <- expomicset |>
#         MultiAssayExperiment::metadata() |>
#         purrr::pluck("enrichment") |>
#         purrr::pluck(geneset)
#
#     if (!is.null(mirna_assays)) {
#         enrich_res <- enrich_res |>
#             dplyr::filter(!exp_name %in% mirna_assays)
#     }
#
#     pc_dfs <- purrr::pmap_dfr(
#         enrich_res |>
#             dplyr::distinct(
#                 exp_name,
#                 Cluster,
#                 go_group
#             ),
#         function(exp_name,
#                  Cluster,
#                  go_group) {
#             df <- enrich_res |>
#                 dplyr::filter(
#                     exp_name == !!exp_name,
#                     Cluster == !!Cluster,
#                     go_group == !!go_group
#                 )
#
#             genes <- unique(unlist(stringr::str_split(df$geneID, "/")))
#
#             if (length(genes) < min_genes) {
#                 return(NULL)
#             }
#
#             exp <- .update_assay_colData(expomicset, exp_name)
#
#             if (!is.null(feature_col)) {
#                 genes <- exp |>
#                     tidybulk::pivot_transcript() |>
#                     dplyr::filter(!!rlang::sym(feature_col) %in% genes) |>
#                     dplyr::pull(.feature)
#             }
#
#             assay <- SummarizedExperiment::assay(exp)
#
#             assay <- assay[rownames(assay) %in% genes, , drop = FALSE]
#
#             if (nrow(assay) < 2 || all(apply(assay, 1, var) == 0)) {
#                 return(NULL)
#             }
#
#             pc1 <- prcomp(t(log2(assay + 1)), scale. = TRUE)$x[, 1]
#
#             id <- rownames(pc1)
#
#             tibble::tibble(
#                 id = id,
#                 !!paste("PC", exp_name, Cluster, go_group, sep = "/") := pc1
#             )
#         }
#     )
#
#     return(pc_dfs)
# }

# --- R2 Function ---------
#' Calculate R-squared and adjusted R-squared from a GLM
#'
#' This internal helper function computes the R-squared and adjusted R-squared
#' values from a fitted generalized linear model (GLM) object.
#'
#' @param model A fitted \code{glm} object.
#' @param n Optional integer. Number of observations used to fit the model.
#'   If \code{NULL}, it is inferred from \code{model$model}.
#'
#' @return A list with two elements:
#'   - \code{r2}: R-squared
#'   - \code{adj_r2}: Adjusted R-squared
#'
#' @details
#' The adjusted R-squared is calculated as:
#' \deqn{1 - \frac{(n - 1)}{(n - p)}(1 - R^2)}
#' where \code{n} is the number of observations and \code{p} is
#' the number of parameters (including the intercept).
#'
#' @keywords internal
#' @noRd
.calc_r2 <- function(
    model,
    model_data,
    family = NULL,
    outcome = NULL) {
    fam <- if (is.null(family)) stats::family(model)$family else family
    n <- nrow(model_data)
    p <- length(stats::coef(model)) # includes intercept

    if (identical(fam, "gaussian")) {
        # regular adjusted R2
        r2 <- 1 - model$deviance / model$null.deviance
        adj_r2 <- 1 - ((n - 1) / (n - p)) * (1 - r2)
        return(list(r2 = r2, adj_r2 = adj_r2))
    }

    if (identical(fam, "binomial")) {
        if (is.null(outcome)) stop("Provide `outcome` for binomial metrics.")
        ll_full <- as.numeric(stats::logLik(model))
        mod_null <- stats::glm(stats::as.formula(paste(outcome, "~ 1")),
            data = model_data, family = stats::binomial()
        )
        ll_null <- as.numeric(stats::logLik(mod_null))

        # McFadden and adjusted McFadden R2
        r2 <- 1 - (ll_full / ll_null)
        adj_r2 <- 1 - ((ll_full - p) / ll_null)

        return(list(r2 = r2, adj_r2 = adj_r2))
    }

    list(r2 = NA_real_, adj_r2 = NA_real_)
}

# .calc_r2 <- function(
#     model,
#     model_data) {
#     # Compute R-squared and adjusted R-squared
#     n <- nrow(model_data)
#     # includes intercept
#     p <- length(coef(model))
#
#     # calculate the r2 and adj r2
#     r2 <- 1 - model$deviance / model$null.deviance
#     adj_r2 <- 1 - ((n - 1) / (n - p)) * (1 - r2)
#
#     return(list(
#         r2 = r2,
#         adj_r2 = adj_r2
#     ))
# }

# --- Annotate the Model Results By Model Source ---------
#' Annotate association results based on feature source
#'
#' This internal helper function adds metadata annotations to the association
#' results table based on the data source (`omics`, `factors`, or `exposures`).
#'
#' For `omics` and `factors`, it:
#' - Extracts experiment names from the `MultiAssayExperiment` and parses them
#'   from the `term` column.
#' - Adds a `category` column representing the experiment name.
#' - For `factors` from DIABLO or RGCCA integration, prefixes the factor names
#'   with the experiment category.
#'
#' For `exposures`, it:
#' - Joins the `codebook` from the `MultiAssayExperiment` metadata to add
#'   variable annotations.
#'
#' For `omics`, it also:
#' - Joins feature-level metadata from `rowData` across all assays.
#'
#' @param results A data frame of association results, with at least a
#'   `term` column.
#' @param source Character string: the feature source (`"omics"`, `"factors"`,
#'   or `"exposures"`).
#' @param expomicset A `MultiAssayExperiment` containing experiments, metadata,
#'   and codebook information.
#'
#' @return A data frame of association results with additional metadata columns
#'   (e.g., `category`, joined feature metadata).
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr mutate case_when select rename inner_join left_join
#' bind_rows
#' @importFrom stringr str_extract
#' @importFrom purrr pluck
#' @importFrom tibble rownames_to_column
#' @importFrom MultiAssayExperiment experiments metadata
#' @importFrom SummarizedExperiment rowData
.annotate_results_by_source <- function(
    results,
    source,
    expomicset) {
    if (source %in% c("omics", "factors")) {
        exp_names <- names(MultiAssayExperiment::experiments(expomicset)) |>
            (\(chr) gsub(" ", "_", chr))()

        matched <- data.frame(
            exp_name = exp_names,
            exp_name_clean = names(MultiAssayExperiment::experiments(expomicset))
        )

        results <- results |>
            dplyr::mutate(
                exp_name = stringr::str_extract(
                    term, paste0("^(", paste(exp_names, collapse = "|"), ")")
                ),
                term = dplyr::case_when(
                    grepl(paste(exp_names, collapse = "|"), term) ~
                        gsub(paste0(
                            "(", paste0(exp_names, "_", collapse = "|"), ")"
                        ), "", term),
                    .default = term
                )
            ) |>
            dplyr::inner_join(matched, by = "exp_name") |>
            dplyr::select(-exp_name) |>
            dplyr::rename(category = exp_name_clean)

        if (source == "factors") {
            method <- MultiAssayExperiment::metadata(expomicset) |>
                purrr::pluck(
                    "multiomics_integration",
                    "integration_results",
                    "method"
                )

            if (method %in% c("DIABLO", "RGCCA")) {
                results <- dplyr::mutate(
                    results,
                    term = paste(category, term, sep = " ")
                )
            }
        }
    }

    if (source == "exposures") {
        results <- results |>
            dplyr::left_join(
                MultiAssayExperiment::metadata(expomicset)$codebook,
                by = c("term" = "variable")
            )
    }

    if (source == "omics") {
        feature_df <- lapply(
            names(MultiAssayExperiment::experiments(expomicset)),
            function(name) {
                SummarizedExperiment::rowData(
                    MultiAssayExperiment::experiments(expomicset)[[name]]
                ) |>
                    as.data.frame() |>
                    tibble::rownames_to_column(".feature") |>
                    dplyr::mutate(exp_name = name)
            }
        ) |>
            dplyr::bind_rows()

        results <- results |>
            dplyr::mutate(category = gsub("_", " ", category)) |>
            dplyr::left_join(
                feature_df,
                by = c("term" = ".feature", "category" = "exp_name")
            )
    }

    return(results)
}
