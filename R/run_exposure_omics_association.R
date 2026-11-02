#' Run Exposure-Omics Association
#'
#' Test associations between each exposure and each omics feature using
#' limma's linear modeling framework.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing exposome
#'   and omics data.
#' @param formula A one-sided formula describing the model to fit for each
#'   exposure. Use `.feature` as a placeholder for the exposure being tested;
#'   it is replaced by each variable in `exposures` in turn. Any other terms
#'   must be columns in `colData` and are treated as covariates. Defaults to
#'   `~ .feature` (no covariates).
#' @param exposures Character vector of exposure variable names to test.
#'   If `NULL`, uses all variables from the codebook.
#' @param exp_name Name(s) of the omics assay(s) to test against. If `NULL`,
#'   uses all assays.
#' @param scaling_method Character. Scaling method to apply before modeling.
#'   Options include `"none"` (default).
#' @param correction_method Method for p-value adjustment. Default is `"fdr"`.
#' @param top_pct Top X% of features to retain using either mean or variance
#'  which is specified by `filter_by`. If `NULL`, no features will be filtered.
#' @param filter_by Determination of how to filter omics features either by
#'   mean or variance.
#' @param action If `"add"` (default), saves results to metadata,
#'   if `"get"`, returns results as a data frame.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`.
#'   Otherwise, returns a tibble with association results.
#'
#' @details
#' This function uses limma to test associations between multiple exposures
#' and omics features. For each exposure, the model given by `formula` is
#' fit with each omics feature as the outcome, e.g.
#'
#' `omics_feature ~ .feature + covariate1 + covariate2 + ...`
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run exposure-omics association
#' mae <- mae |>
#'     run_exposure_omics_association(
#'         formula = ~ .feature + age + sex,
#'         exposures = c("exposure_pm25", "exposure_no2")
#'     )
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix as.formula
#' @export
run_exposure_omics_association <- function(
  exposomicset,
  formula = ~.feature,
  exposures = NULL,
  exp_name = NULL,
  scaling_method = "none",
  correction_method = "fdr",
  top_pct = NULL,
  filter_by = c("variance", "mean"),
  action = "add"
) {
    .check_suggested("limma")
    filter_by <- match.arg(filter_by)

    # Validate formula
    if (!inherits(formula, "formula") || length(formula) != 2) {
        stop("`formula` must be a one-sided formula, e.g. ~ .feature + age")
    }
    if (!".feature" %in% all.vars(formula)) {
        stop("`formula` must include the `.feature` placeholder, e.g. ~ .feature + age")
    }

    # Covariates are every term in the formula other than the placeholder
    covariates <- setdiff(all.vars(formula), ".feature")

    # Get exposure names from codebook if not provided
    if (is.null(exposures)) {
        exposures <- MultiAssayExperiment::metadata(exposomicset)$codebook$variable
    }

    # Validate exposures and covariates exist
    missing_vars <- exposomicset |>
        MultiAssayExperiment::colData() |>
        colnames() |>
        (\(x) setdiff(c(exposures, covariates), x))()

    if (length(missing_vars) > 0) {
        stop("Variables not found in colData: ", paste(missing_vars, collapse = ", "))
    }

    # Remove covariates from exposures if present
    exposures <- setdiff(exposures, covariates)

    # Get assay names
    if (is.null(exp_name)) {
        exp_name <- names(MultiAssayExperiment::experiments(exposomicset))
    }

    message(sprintf(
        "Testing %d exposures across %d assays using model: %s",
        length(exposures), length(exp_name), deparse(formula)
    ))

    # Run analysis across all assays
    results <- exp_name |>
        lapply(function(assay_name) {
            message("Processing assay: ", assay_name)

            exposomicset |>
                .update_assay_colData(assay_name) |>
                .filter_top_features(
                    top_pct = top_pct,
                    filter_by = filter_by
                ) |>
                .run_exposure_omics_limma(
                    exposures = exposures,
                    formula = formula,
                    scaling_method = scaling_method
                )
        }) |>
        stats::setNames(exp_name) |>
        dplyr::bind_rows(.id = "exp_name") |>
        dplyr::filter(contrast == ".feature") |>
        dplyr::mutate(p_adjust = stats::p.adjust(p.value, method = correction_method)) |>
        dplyr::select(-contrast) |>
        tibble::as_tibble()

    # Add exposure annotations from codebook
    codebook <- MultiAssayExperiment::metadata(exposomicset)$codebook
    if (!is.null(codebook) && "variable" %in% colnames(codebook)) {
        results <- results |>
            dplyr::left_join(codebook, by = c("exposure" = "variable"))
    }

    if (action == "add") {
        MultiAssayExperiment::metadata(exposomicset)$association$exposure_omics <- list(
            results_df = results,
            formula = formula,
            covariates = covariates,
            exposures = exposures,
            exp_name = exp_name,
            method = "limma_trend"
        )

        step_record <- list(
            run_exposure_omics_association = list(
                timestamp = Sys.time(),
                params = list(
                    formula = deparse(formula),
                    exposures = exposures,
                    exp_name = exp_name,
                    scaling_method = scaling_method,
                    correction_method = correction_method
                ),
                notes = sprintf(
                    "Tested %d exposures against %d assays using limma-trend (%s)",
                    length(exposures), length(exp_name), deparse(formula)
                )
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else {
        return(results)
    }
}

#' Run limma-trend for exposure-omics associations
#' @keywords internal
#' @noRd
.run_exposure_omics_limma <- function(
  se,
  exposures,
  formula,
  scaling_method
) {
    exposures |>
        lapply(function(exp_var) {
            # Scale numeric exposures, keep categorical as is
            if (is.numeric(SummarizedExperiment::colData(se)[[exp_var]])) {
                SummarizedExperiment::colData(se)$.feature <- as.numeric(
                    scale(SummarizedExperiment::colData(se)[[exp_var]])
                )
            } else {
                SummarizedExperiment::colData(se)$.feature <-
                    SummarizedExperiment::colData(se)[[exp_var]]
            }

            .run_limma_trend(
                se = se,
                formula = formula,
                abundance_col = SummarizedExperiment::assayNames(se)[1],
                scaling_method = scaling_method
            ) |>
                dplyr::mutate(exposure = exp_var) |>
                dplyr::rename(
                    feature = feature,
                    estimate = logFC,
                    t_statistic = t,
                    p.value = P.Value,
                    p_adjust = adj.P.Val
                )
        }) |>
        dplyr::bind_rows()
}

#' Filter to top features by variance or mean
#' @keywords internal
#' @noRd
.filter_top_features <- function(
  se,
  top_pct = NULL,
  filter_by = "variance"
) {
    if (is.null(top_pct)) {
        return(se)
    }

    mat <- SummarizedExperiment::assay(se)

    scores <- switch(filter_by,
        variance = apply(mat, 1, var, na.rm = TRUE),
        mean = apply(mat, 1, mean, na.rm = TRUE)
    )

    n_features <- ceiling(nrow(se) * top_pct / 100)
    top_features <- names(sort(scores, decreasing = TRUE))[seq_len(n_features)]

    message(sprintf(
        "Filtered to top %d%% (%d features) by %s",
        top_pct, n_features, filter_by
    ))

    se[top_features, ]
}
