#' Run Association Analysis
#'
#' Perform GLM-based association testing between a specified outcome and
#' features from exposures, omics, or latent factors. Each feature is tested
#' in its own model, built from a user-supplied formula. Models can be run on
#' subsets of samples, stratified by a grouping variable, or repeated across
#' several formulas. Every fitted model is stored under a name so that
#' specifications can be compared later with `compare_association()`.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing data
#' and metadata.
#' @param formula A two-sided formula whose left-hand side is a single
#' column of `colData` (the outcome) and whose right-hand side contains the
#' covariates, e.g. `smoker ~ age * sex`. Each feature is added to the model
#' as a main effect. To control how the feature enters the model (for example
#' to test an interaction), use the placeholder `.feature` on the right-hand
#' side, e.g. `smoker ~ .feature * sex + age`. Use `outcome ~ 1` for an
#' unadjusted model. Alternatively, a **named list of formulas**
#' (e.g. `list(crude = y ~ 1, adjusted = y ~ age + sex)`) fits every
#' specification in one call, each stored under its list name.
#' @param source Source of features to test. One of `"omics"`,
#' `"exposures"`, `"factors"`.
#' @param feature_set Optional character vector of exposure or GO terms to test.
#' @param log_trans Optional boolean value dictating whether or not to log2
#' transform omics features.
#' @param top_n Optional integer: if using omics source, select top `n`
#'  most variable features.
#' @param family GLM family; `"gaussian"` or `"binomial"`.
#' @param correction_method Method for p-value adjustment (default: `"fdr"`).
#' Adjustment is performed separately for each model and each effect
#' (see Value).
#' @param subset Optional bare expression evaluated against the (unscaled)
#' `colData` that selects the samples to model, e.g. `sex == "female"` or
#' `age > 40 & smoker == "no"`.
#' @param strata Optional single column name of `colData`. The model is fit
#' separately within each level and stored as `[MODEL_NAME]_[LEVEL_IN_STRATA]`. The stratifying
#' variable must not appear in `formula`.
#' @param name Optional name for the stored model. Defaults to `"default"`,
#' or to the text of `subset` when one is given. Ignored when `formula` is a
#' named list (the list names are used).
#' @param overwrite Logical; if `TRUE`, replace an existing model with the same
#' name instead of erroring or appending features to it (default: `FALSE`).
#' @param action If `"add"` (default), saves results to metadata;
#' else returns results as list.
#' @param feature_col The column in `rowData`
#' for matching gene symbols or IDs.
#' @param mirna_assays Optional character vector of assays to exclude
#' when extracting GO terms.
#'
#' @details
#' Numeric `colData` columns are z-scored once on the **full** sample before any
#' subsetting, so effect sizes are on the same standardized scale in every
#' subset and stratum and can be compared directly.
#'
#' If a model with the same `name` already exists for the source and was fit
#' with the same formula, subset, and family, new features are appended to it
#' and p-values are re-adjusted across the combined set. A different
#' specification under an existing name is an error unless `overwrite = TRUE`.
#'
#' With `strata`, a stratum that cannot be fit (for example, a single outcome
#' level or a constant covariate) is skipped with a warning.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`.
#'  Otherwise, returns a list of:
#' - `models`: a named list of fitted models, each containing `results_df`
#'   (tidy summary of associations), `formula`, `covariates`, `family`,
#'   `subset`, `strata_var`, `stratum`, and `n_samples`
#' - `model_data`: the (full) data used in the GLMs
#'
#' In `results_df`, `term` is the feature, `coef` is the model coefficient
#' name, `effect` is `"main"` for the feature's main effect or the name of the
#' interacting term (e.g. `"sexMale"`), and `n` / `n_events` give the
#' observations (and binomial events) used.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # covariate-adjusted main effect
#' mae <- mae |>
#'     run_association(
#'         formula = smoker ~ age,
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial"
#'     )
#'
#' # several specifications at once
#' mae <- mae |>
#'     run_association(
#'         formula = list(
#'             crude = smoker ~ 1,
#'             adjusted = smoker ~ age
#'         ),
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial"
#'     )
#'
#' \dontrun{
#' # subset of samples (assumes a `sex` column in colData)
#' mae <- mae |>
#'     run_association(
#'         formula = smoker ~ age,
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial",
#'         subset = sex == "female",
#'         name = "female"
#'     )
#'
#' # stratified by sex, then compared
#' mae <- mae |>
#'     run_association(
#'         formula = smoker ~ age,
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial",
#'         strata = "sex",
#'         name = "by_sex"
#'     )
#'
#' mae |>
#'     compare_association(
#'         source = "exposures",
#'         models = c("by_sex_female", "by_sex_male")
#'     )
#' }
#' @importFrom stats binomial glm logLik as.formula coef family p.adjust
#' @importFrom rlang .data
#' @export
run_association <- function(
  exposomicset,
  formula,
  source = c("omics", "exposures", "factors"),
  feature_set = NULL,
  log_trans = TRUE,
  top_n = NULL,
  family = "gaussian",
  correction_method = "fdr",
  subset = NULL,
  strata = NULL,
  name = NULL,
  overwrite = FALSE,
  action = "add",
  feature_col = NULL,
  mirna_assays = NULL
) {
    .check_suggested(pkg = "tidyselect")
    source <- match.arg(source)

    subset_quo <- rlang::enquo(subset)
    subset_txt <- if (rlang::quo_is_null(subset_quo)) {
        NULL
    } else {
        rlang::quo_text(subset_quo)
    }

    formulas <- .as_formula_list(formula, name, subset_txt)
    parsed_all <- purrr::map(formulas, .parse_association_formula)
    outcomes <- unique(purrr::map_chr(parsed_all, "outcome"))

    if (!is.null(strata)) {
        purrr::walk(parsed_all, function(p) {
            if (strata %in% p$covariates) {
                stop(sprintf(
                    "Stratifying variable '%s' also appears in `formula`; remove it from the formula.",
                    strata
                ))
            }
        })
    }

    raw <- exposomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame()

    # scale numeric columns once on the full sample (comparable across subsets)
    data <- raw |>
        dplyr::mutate(dplyr::across(
            tidyselect::where(is.numeric) & !tidyselect::all_of(outcomes),
            ~ as.numeric(scale(.))
        ))

    sample_sets <- .build_sample_sets(raw, subset_quo, subset_txt, strata)

    features_df <- switch(source,
        omics = .extract_omics_features(
            exposomicset,
            log_trans = log_trans,
            top_n
        ),
        exposures = .extract_exposures(
            data,
            feature_set
        ),
        factors = .extract_latent_factors(exposomicset)
    )

    # create the model data
    model_data <- if (source %in% c("exposures", "factors")) {
        data
    } else {
        dplyr::left_join(
            tibble::rownames_to_column(data, "id"),
            features_df,
            by = "id"
        ) |>
            tibble::column_to_rownames("id")
    }

    # one spec per formula x sample set
    specs <- purrr::imap(formulas, function(fml, base_name) {
        purrr::map(sample_sets, function(ss) {
            list(
                model_name = if (is.null(ss$level)) {
                    base_name
                } else {
                    paste(base_name, ss$level, sep = "_")
                },
                base_name = base_name,
                formula = fml,
                ss = ss
            )
        })
    }) |>
        unlist(recursive = FALSE, use.names = FALSE)

    message("Running GLMs.")

    models <- specs |>
        purrr::map(function(s) {
            tryCatch(
                .fit_association_model(
                    spec = s,
                    model_data = model_data,
                    feature_cols = setdiff(colnames(features_df), "id"),
                    family = family,
                    correction_method = correction_method,
                    source = source,
                    exposomicset = exposomicset
                ),
                error = function(e) {
                    if (is.null(strata)) {
                        stop(conditionMessage(e), call. = FALSE)
                    }
                    warning(
                        sprintf("Skipping '%s': %s", s$model_name, conditionMessage(e)),
                        call. = FALSE
                    )
                    NULL
                }
            )
        }) |>
        purrr::set_names(purrr::map_chr(specs, "model_name")) |>
        purrr::compact()

    if (length(models) == 0) {
        stop("No models could be fit.")
    }

    if (action == "add") {
        key <- paste0("assoc_", source)

        store <- MultiAssayExperiment::metadata(exposomicset)$association[[key]] |>
            .normalize_assoc_store()

        for (nm in names(models)) {
            store$models[[nm]] <- .merge_association_model(
                old = store$models[[nm]],
                new = models[[nm]],
                name = nm,
                source = source,
                overwrite = overwrite
            )
        }

        MultiAssayExperiment::metadata(exposomicset)$association[[key]] <- store

        # Add in a step record
        step_record <- list(
            run_association = list(
                timestamp = Sys.time(),
                params = list(
                    formula = purrr::map_chr(formulas, .formula_key),
                    models = names(models),
                    source = source,
                    feature_set = feature_set,
                    top_n = top_n,
                    family = family,
                    correction_method = correction_method,
                    subset = subset_txt,
                    strata = strata,
                    feature_col = feature_col,
                    mirna_assays = mirna_assays
                ),
                notes = paste0(
                    "Performed association analysis using source: ",
                    source
                )
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    }

    list(
        models = models,
        model_data = model_data
    )
}


# --- Formula Helpers ----------
#' Parse and validate an association formula
#'
#' @param formula A two-sided formula.
#' @return A list with `outcome`, `covariates`, and `has_placeholder`.
#' @keywords internal
#' @noRd
.parse_association_formula <- function(formula) {
    if (!inherits(formula, "formula") || length(formula) != 3) {
        stop("`formula` must be two-sided, e.g. `outcome ~ age * sex`.")
    }
    if (!is.name(formula[[2]])) {
        stop("The left-hand side of `formula` must be a single column name.")
    }

    rhs_vars <- all.vars(formula[[3]])
    if ("." %in% rhs_vars) {
        stop("`.` is not supported in `formula`; list covariates explicitly.")
    }

    list(
        outcome = as.character(formula[[2]]),
        covariates = setdiff(rhs_vars, ".feature"),
        has_placeholder = ".feature" %in% rhs_vars
    )
}

#' Insert a feature into the user formula
#'
#' If the formula contains `.feature`, it is replaced by the feature name.
#' Otherwise the feature is appended as a main effect.
#'
#' @keywords internal
#' @noRd
.build_feature_formula <- function(formula, fcol, has_placeholder) {
    feature <- rlang::sym(fcol)
    rhs <- formula[[3]]

    rhs <- if (has_placeholder) {
        do.call(substitute, list(rhs, list(.feature = feature)))
    } else {
        rlang::expr(!!rhs + !!feature)
    }

    rlang::new_formula(formula[[2]], rhs, env = environment(formula))
}

#' Single-string key for comparing / recording formulas
#' @keywords internal
#' @noRd
.formula_key <- function(formula) {
    paste(deparse(formula), collapse = "")
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
.extract_omics_features <- function(exposomicset,
                                    log_trans = TRUE,
                                    top_n = NULL) {
    if (log_trans) {
        log2_assays <- .log2_multiassay(exposomicset)
    } else {
        log2_assays <- exposomicset
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
# moved to run_multi_omics_integration
# --- R2 Function ---------
#' Calculate R-squared and adjusted R-squared from a GLM
#'
#' This internal helper function computes the R-squared and adjusted R-squared
#' values from a fitted generalized linear model (GLM) object.
#'
#' @param model A fitted \code{glm} object.
#' @param model_data The data frame used to fit the model.
#' @param family Optional GLM family name.
#' @param outcome Outcome variable name (required for binomial metrics).
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
  outcome = NULL
) {
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
#' @param exposomicset A `MultiAssayExperiment` containing experiments, metadata,
#'   and codebook information.
#'
#' @return A data frame of association results with additional metadata columns
#'   (e.g., `category`, joined feature metadata).
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr mutate case_when select rename inner_join left_join bind_rows
#' @importFrom stringr str_extract
#' @importFrom purrr pluck
#' @importFrom tibble rownames_to_column
#' @importFrom MultiAssayExperiment experiments metadata
#' @importFrom SummarizedExperiment rowData
.annotate_results_by_source <- function(
  results,
  source,
  exposomicset
) {
    if (source == "factors") {
        method <- MultiAssayExperiment::metadata(exposomicset) |>
            purrr::pluck("multiomics_integration", "integration_results", "method")

        if (method %in% c("DIABLO", "RGCCA")) {
            # These have assay prefixes in term names
            exp_names <- names(MultiAssayExperiment::experiments(exposomicset)) |>
                (\(chr) gsub(" ", "_", chr))()
            matched <- data.frame(
                exp_name = exp_names,
                exp_name_clean = names(MultiAssayExperiment::experiments(exposomicset))
            )
            results <- results |>
                dplyr::mutate(
                    exp_name = stringr::str_extract(
                        term, paste0("^(", paste(exp_names, collapse = "|"), ")")
                    ),
                    term = dplyr::case_when(
                        grepl(paste(exp_names, collapse = "|"), term) ~
                            gsub(paste0("(", paste0(exp_names, "_", collapse = "|"), ")"), "", term),
                        .default = term
                    )
                ) |>
                dplyr::inner_join(matched, by = "exp_name") |>
                dplyr::select(-exp_name) |>
                dplyr::rename(category = exp_name_clean)

            results <- dplyr::mutate(results, term = paste(category, term, sep = " "))
        }
        # MOFA, MCIA, MCCA have no assay prefix, so no exp_name/category

        return(results)
    }

    if (source == "omics") {
        exp_names <- names(MultiAssayExperiment::experiments(exposomicset)) |>
            (\(chr) gsub(" ", "_", chr))()
        matched <- data.frame(
            exp_name = exp_names,
            exp_name_clean = names(MultiAssayExperiment::experiments(exposomicset))
        )
        results <- results |>
            dplyr::mutate(
                exp_name = stringr::str_extract(
                    term, paste0("^(", paste(exp_names, collapse = "|"), ")")
                ),
                term = dplyr::case_when(
                    grepl(paste(exp_names, collapse = "|"), term) ~
                        gsub(paste0("(", paste0(exp_names, "_", collapse = "|"), ")"), "", term),
                    .default = term
                )
            ) |>
            dplyr::inner_join(matched, by = "exp_name") |>
            dplyr::select(-exp_name) |>
            dplyr::rename(category = exp_name_clean) |>
            dplyr::mutate(category = gsub("_", " ", category))

        # Add feature metadata
        feature_df <- lapply(
            names(MultiAssayExperiment::experiments(exposomicset)),
            function(name) {
                SummarizedExperiment::rowData(
                    MultiAssayExperiment::experiments(exposomicset)[[name]]
                ) |>
                    as.data.frame() |>
                    tibble::rownames_to_column(".feature") |>
                    dplyr::mutate(exp_name = name)
            }
        ) |>
            dplyr::bind_rows()

        results <- results |>
            dplyr::left_join(feature_df, by = c("term" = ".feature", "category" = "exp_name"))

        return(results)
    }

    if (source == "exposures") {
        results <- results |>
            dplyr::left_join(
                MultiAssayExperiment::metadata(exposomicset)$codebook,
                by = c("term" = "variable")
            )
        return(results)
    }

    results
}
