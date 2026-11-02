# --- Model Store Helpers ----------
# Results live at metadata(mae)$association$assoc_<source>$models[[NAME_OF_MODEL]].
# Older objects stored a single model directly in assoc_<source>; those are
# transparently treated as a model called "default".

#' Normalise a stored association object into `list(models = ...)`
#' @keywords internal
#' @noRd
.normalize_assoc_store <- function(store) {
    if (is.null(store)) {
        return(list(models = list()))
    }

    # legacy: a single model stored directly
    if (is.null(store$models)) {
        store <- list(models = list(default = store))
    }

    store$models <- purrr::map(store$models, function(m) {
        if (!"effect" %in% names(m$results_df)) {
            m$results_df$effect <- "main"
        }
        if (is.null(m$formula) && !is.null(m$results_df$outcome)) {
            m$formula <- stats::reformulate(
                if (length(m$covariates) > 0) m$covariates else "1",
                response = m$results_df$outcome[[1]]
            )
        }
        m
    })

    store
}

#' Fetch one stored model
#' @keywords internal
#' @noRd
.get_assoc_model <- function(exposomicset, source, model = NULL, arg = "model") {
    store <- MultiAssayExperiment::metadata(exposomicset)$association[[
        paste0("assoc_", source)
    ]] |>
        .normalize_assoc_store()

    if (length(store$models) == 0) {
        stop(
            "Association results for source '",
            source,
            "' not found in metadata. Run `run_association()` first."
        )
    }

    if (is.null(model)) {
        if (length(store$models) > 1) {
            stop(
                "Multiple models stored for source '", source, "': ",
                paste(names(store$models), collapse = ", "),
                ". Specify `", arg, "`."
            )
        }
        return(store$models[[1]])
    }

    if (!model %in% names(store$models)) {
        stop(
            "Model '", model, "' not found for source '", source, "'. Available: ",
            paste(names(store$models), collapse = ", ")
        )
    }

    store$models[[model]]
}

#' Fetch stored results for a source, or NULL when none exist
#'
#' Used by plots that combine several sources (e.g. `plot_manhattan()`), where
#' a missing source is allowed.
#' @keywords internal
#' @noRd
.get_assoc_results <- function(exposomicset, source, model = NULL, arg = "model") {
    store <- MultiAssayExperiment::metadata(exposomicset)$association[[
        paste0("assoc_", source)
    ]] |>
        .normalize_assoc_store()

    if (length(store$models) == 0) {
        return(NULL)
    }

    .get_assoc_model(exposomicset, source, model, arg)$results_df
}

#' Combine a newly fit model with a stored one of the same name
#' @keywords internal
#' @noRd
.merge_association_model <- function(old, new, name, source, overwrite = FALSE) {
    if (is.null(old) || overwrite) {
        return(new)
    }

    same_spec <- identical(.formula_key(old$formula), .formula_key(new$formula)) &&
        identical(old$subset, new$subset) &&
        (is.null(old$family) || identical(old$family, new$family))

    if (!same_spec) {
        stop(
            "Model '", name, "' already exists for source '", source,
            "' with a different specification. Use a different `name` ",
            "or set `overwrite = TRUE`."
        )
    }

    if (any(new$results_df$term %in% old$results_df$term)) {
        stop(
            "Association results for this feature already exist in model '",
            name, "'. Set `overwrite = TRUE` to replace them."
        )
    }

    new$results_df <- dplyr::bind_rows(old$results_df, new$results_df) |>
        .adjust_pvalues(new$correction_method)

    new
}

#' Adjust p-values within each effect
#' @keywords internal
#' @noRd
.adjust_pvalues <- function(df, method) {
    df |>
        dplyr::group_by(effect) |>
        dplyr::mutate(p_adjust = stats::p.adjust(p.value, method = method)) |>
        dplyr::ungroup()
}

#' Name the formula(s) passed to run_association()
#' @keywords internal
#' @noRd
.as_formula_list <- function(formula, name, subset_txt) {
    if (inherits(formula, "formula")) {
        base <- if (!is.null(name)) {
            name
        } else if (!is.null(subset_txt)) {
            subset_txt
        } else {
            "default"
        }
        return(stats::setNames(list(formula), base))
    }

    valid_list <- is.list(formula) &&
        length(formula) > 0 &&
        all(vapply(formula, inherits, logical(1), what = "formula")) &&
        !is.null(names(formula)) &&
        all(nzchar(names(formula))) &&
        !anyDuplicated(names(formula))

    if (!valid_list) {
        stop("`formula` must be a formula or a uniquely named list of formulas.")
    }

    if (!is.null(name)) {
        warning("`name` is ignored when `formula` is a named list.", call. = FALSE)
    }

    formula
}

#' Build the sample sets (all / subset / strata) to model
#' @keywords internal
#' @noRd
.build_sample_sets <- function(raw, subset_quo, subset_txt, strata) {
    ids <- rownames(raw)

    if (!rlang::quo_is_null(subset_quo)) {
        keep <- rlang::eval_tidy(subset_quo, data = raw)

        if (!is.logical(keep) || length(keep) != nrow(raw)) {
            stop("`subset` must evaluate to one logical value per sample.")
        }

        ids <- ids[which(keep)]

        if (length(ids) == 0) {
            stop("`subset` matched no samples.")
        }
    }

    if (is.null(strata)) {
        return(list(list(
            level = NULL,
            ids = ids,
            subset = subset_txt,
            strata_var = NULL
        )))
    }

    if (length(strata) != 1 || !strata %in% names(raw)) {
        stop("`strata` must be a single column name of colData.")
    }

    vals <- raw[ids, strata]
    levs <- sort(unique(stats::na.omit(vals)))

    lapply(levs, function(lv) {
        list(
            level = as.character(lv),
            ids = ids[which(vals == lv)],
            subset = subset_txt,
            strata_var = strata
        )
    })
}

#' Label for the samples a model was fit on ("all", the subset, and/or stratum)
#' @keywords internal
#' @noRd
.sample_key <- function(ss) {
    parts <- c(
        ss$subset,
        if (!is.null(ss$level)) paste0(ss$strata_var, "=", ss$level)
    )

    if (length(parts) == 0) "all" else paste(parts, collapse = " & ")
}

#' Fit all features for one model spec
#' @keywords internal
#' @noRd
.fit_association_model <- function(
  spec,
  model_data,
  feature_cols,
  family,
  correction_method,
  source,
  exposomicset
) {
    parsed <- .parse_association_formula(spec$formula)
    outcome <- parsed$outcome
    covariates <- parsed$covariates

    missing_vars <- setdiff(c(outcome, covariates), colnames(model_data))
    if (length(missing_vars) > 0) {
        stop(
            "Variables in `formula` not found in colData: ",
            paste(missing_vars, collapse = ", ")
        )
    }

    model_data <- model_data[
        rownames(model_data) %in% spec$ss$ids, ,
        drop = FALSE
    ]

    if (nrow(model_data) == 0) {
        stop("No samples available for this model.")
    }

    # ensure outcome variable is factor if binomial
    if (family == "binomial") {
        model_data[[outcome]] <- as.factor(model_data[[outcome]])

        if (length(unique(stats::na.omit(model_data[[outcome]]))) != 2) {
            stop("Binary outcome must have exactly 2 levels")
        }
    }

    # ensure covariates have more than one value
    purrr::walk(covariates, function(cov) {
        if (length(unique(stats::na.omit(model_data[[cov]]))) <= 1) {
            stop(sprintf(
                "Covariate %s must have more than one unique value.",
                cov
            ))
        }
    })

    results <- purrr::map_dfr(feature_cols, function(fcol) {
        model <- tryCatch(
            glm(
                .build_feature_formula(
                    spec$formula,
                    fcol,
                    parsed$has_placeholder
                ),
                data = model_data,
                family = family
            ),
            error = function(e) NULL
        )

        if (is.null(model)) {
            return(NULL)
        }

        fit_r2 <- .calc_r2(
            model,
            model_data,
            family = family,
            outcome = outcome
        )

        # keep the feature's main effect and any interactions involving it
        broom::tidy(model) |>
            dplyr::filter(purrr::map_lgl(
                strsplit(term, ":", fixed = TRUE),
                ~ fcol %in% .x
            )) |>
            dplyr::mutate(
                coef = term,
                effect = purrr::map_chr(
                    strsplit(coef, ":", fixed = TRUE),
                    ~ paste(setdiff(.x, fcol), collapse = ":")
                ),
                effect = dplyr::if_else(effect == "", "main", effect),
                term = fcol,
                r2 = fit_r2$r2,
                adj_r2 = fit_r2$adj_r2,
                n = stats::nobs(model),
                n_events = if (family == "binomial") sum(model$y) else NA_real_
            )
    })

    if (nrow(results) == 0) {
        stop("No models could be fit.")
    }

    results <- results |>
        dplyr::mutate(
            outcome = outcome,
            model = spec$model_name
        ) |>
        .adjust_pvalues(correction_method) |>
        .annotate_results_by_source(
            source = source,
            exposomicset = exposomicset
        )

    list(
        results_df = results,
        formula = spec$formula,
        covariates = covariates,
        family = family,
        subset = spec$ss$subset,
        strata_var = spec$ss$strata_var,
        stratum = spec$ss$level,
        base_name = spec$base_name,
        sample_key = .sample_key(spec$ss),
        n_samples = nrow(model_data),
        correction_method = correction_method
    )
}


# --- Model Listing ----------
#' List Stored Association Models
#'
#' Summarise every association model stored in a `MultiAssayExperiment`,
#' across all sources.
#'
#' @param exposomicset A `MultiAssayExperiment` with results from
#' [run_association()].
#'
#' @return A tibble with one row per stored model: `source`, `model`,
#' `formula`, `subset`, `stratum`, `family`, `n_samples`, and `n_features`.
#'
#' @examples
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' ) |>
#'     run_association(
#'         formula = smoker ~ age,
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial"
#'     )
#'
#' association_models(mae)
#' @export
association_models <- function(exposomicset) {
    assoc <- MultiAssayExperiment::metadata(exposomicset)$association

    if (is.null(assoc)) {
        return(tibble::tibble())
    }

    or_na <- function(x) if (is.null(x)) NA_character_ else as.character(x)

    purrr::imap_dfr(assoc, function(obj, key) {
        purrr::imap_dfr(
            .normalize_assoc_store(obj)$models,
            function(m, nm) {
                tibble::tibble(
                    source = sub("^assoc_", "", key),
                    model = nm,
                    formula = or_na(if (!is.null(m$formula)) .formula_key(m$formula)),
                    subset = or_na(m$subset),
                    stratum = or_na(m$stratum),
                    family = or_na(m$family),
                    n_samples = if (is.null(m$n_samples)) NA_integer_ else m$n_samples,
                    n_features = dplyr::n_distinct(m$results_df$term)
                )
            }
        )
    })
}


# --- Model Comparison ----------
#' Compare Association Models
#'
#' Line up the results of two or more stored models (e.g. crude vs adjusted,
#' or female vs male strata) feature by feature.
#'
#' @param exposomicset A `MultiAssayExperiment` with results from
#' [run_association()].
#' @param source Association source: `"omics"`, `"exposures"`, `"factors"`,
#' or `"go_pcs"`.
#' @param models Character vector of stored model names to compare. Default
#' `NULL` compares every model stored for the source.
#' @param reference Model to use as the reference for the derived columns.
#' Defaults to the first of `models`.
#' @param effect Which effect to compare: `"main"` (default) or an interacting
#' coefficient name such as `"sexMale"`.
#' @param terms Optional character vector of features to keep.
#'
#' @details
#' The `relation` column says how each model relates to the reference, which
#' determines which derived statistics are filled in:
#' - `"same_sample"`: same samples, different model (e.g. crude vs adjusted).
#'   `pct_change` is the change-in-estimate,
#'   `100 * (estimate - reference estimate) / reference estimate`;
#'   negative values mean the effect shrank in magnitude.
#' - `"sibling_strata"`: different strata from the same `strata` run.
#'   `z_diff` and `p_diff` test the difference between the two estimates,
#'   `(b1 - b2) / sqrt(se1^2 + se2^2)`, which assumes the strata are
#'   independent. Comparisons are pairwise against the reference.
#' - `"unrelated"`: no derived statistics.
#'
#' @return A tibble with one row per feature and model: `term`, `model`,
#' `relation`, `estimate`, `std.error`, `p.value`, any of `p_adjust`, `n`,
#' `n_events`, then `ref_estimate`, `ref_se`, `reference`, `pct_change`,
#' `z_diff`, and `p_diff`.
#'
#' @examples
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' ) |>
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
#' compare_association(
#'     mae,
#'     source = "exposures"
#' )
#' @export
compare_association <- function(
  exposomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  models = NULL,
  reference = NULL,
  effect = "main",
  terms = NULL
) {
    source <- match.arg(source)

    store <- MultiAssayExperiment::metadata(exposomicset)$association[[
        paste0("assoc_", source)
    ]] |>
        .normalize_assoc_store()

    models <- if (is.null(models)) names(store$models) else models

    unknown <- setdiff(models, names(store$models))
    if (length(unknown) > 0) {
        stop(
            "Model(s) not found: ", paste(unknown, collapse = ", "),
            ". See `association_models()`."
        )
    }
    if (length(models) < 2) {
        stop("Need at least two models to compare. See `association_models()`.")
    }

    reference <- if (is.null(reference)) models[[1]] else reference
    if (!reference %in% models) {
        stop("`reference` must be one of `models`.")
    }

    ref_meta <- store$models[[reference]]

    relations <- purrr::map_chr(models, function(nm) {
        m <- store$models[[nm]]

        if (nm == reference) {
            "reference"
        } else if (identical(ref_meta$sample_key, m$sample_key)) {
            "same_sample"
        } else if (
            !is.null(ref_meta$strata_var) &&
                identical(ref_meta$strata_var, m$strata_var) &&
                identical(ref_meta$base_name, m$base_name) &&
                identical(ref_meta$subset, m$subset)
        ) {
            "sibling_strata"
        } else {
            "unrelated"
        }
    }) |>
        purrr::set_names(models)

    long <- purrr::map_dfr(models, function(nm) {
        store$models[[nm]]$results_df |>
            dplyr::filter(effect == !!effect) |>
            dplyr::mutate(model = nm) |>
            dplyr::select(
                model,
                term,
                estimate,
                std.error,
                p.value,
                dplyr::any_of(c("p_adjust", "n", "n_events"))
            )
    }) |>
        dplyr::filter(is.null(terms) | term %in% terms)

    ref_est <- long |>
        dplyr::filter(model == reference) |>
        dplyr::select(
            term,
            ref_estimate = estimate,
            ref_se = std.error
        )

    long |>
        dplyr::left_join(ref_est, by = "term") |>
        dplyr::mutate(
            reference = reference,
            relation = unname(relations[model]),
            pct_change = dplyr::if_else(
                relation == "same_sample",
                100 * (estimate - ref_estimate) / ref_estimate,
                NA_real_
            ),
            z_diff = dplyr::if_else(
                relation == "sibling_strata",
                (estimate - ref_estimate) / sqrt(std.error^2 + ref_se^2),
                NA_real_
            ),
            p_diff = 2 * stats::pnorm(-abs(z_diff)),
            model = factor(model, levels = models)
        ) |>
        dplyr::arrange(term, model) |>
        dplyr::select(term, model, relation, dplyr::everything())
}
