#' Run Exposure-Omics Association
#'
#' Test associations between each exposure and each omics feature using
#' limma's linear modeling framework.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing exposome
#'   and omics data.
#' @param exposures Character vector of exposure variable names to test.
#'   If `NULL`, uses all variables from the codebook.
#' @param omics_assay Name(s) of the omics assay(s) to test against. If `NULL`,
#'   uses all assays.
#' @param covariates Optional character vector of covariate names to include
#'   in the model.
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
#' and omics features. For each exposure, a linear model is fit with the
#' exposure as the predictor and each omics feature as the outcome,
#' adjusting for covariates.
#'
#' `omics_feature ~ exposure + covariate1 + covariate2 + ...`
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
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         covariates = c("age", "sex")
#'     )
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix as.formula
#' @export
run_exposure_omics_association <- function(
  exposomicset,
  exposures = NULL,
  omics_assay = NULL,
  covariates = NULL,
  scaling_method = "none",
  correction_method = "fdr",
  top_pct = NULL,
  filter_by = c("variance", "mean"),
  action = "add"
) {
    .check_suggested("limma")
    filter_by <- match.arg(filter_by)

    # Get exposure names from codebook if not provided
    if (is.null(exposures)) {
        exposures <- MultiAssayExperiment::metadata(exposomicset)$codebook$variable
    }

    # Validate exposures exist
    col_data <- exposomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame()

    missing_exp <- setdiff(exposures, colnames(col_data))
    if (length(missing_exp) > 0) {
        stop("Exposures not found in colData: ", paste(missing_exp, collapse = ", "))
    }

    # Remove covariates from exposures if present
    exposures <- setdiff(exposures, covariates)

    # Get assay names
    if (is.null(omics_assay)) {
        omics_assay <- names(MultiAssayExperiment::experiments(exposomicset))
    }

    message(sprintf(
        "Testing %d exposures across %d assays",
        length(exposures), length(omics_assay)
    ))

    # Run analysis across all assays
    results <- omics_assay |>
        lapply(function(assay_name) {
            message("Processing assay: ", assay_name)

            exposomicset |>
                .update_assay_colData(assay_name) |>
                .filter_top_features(
                    se,
                    top_n = top_n,
                    top_pct = top_pct,
                    filter_by = filter_by
                ) |>
                .run_exposure_omics_limma(
                    exposures = exposures,
                    covariates = covariates,
                    scaling_method = scaling_method
                )
        }) |>
        stats::setNames(omics_assay) |>
        dplyr::bind_rows(.id = "omics_assay") |>
        dplyr::filter(contrast == ".exposure") |>
        dplyr::mutate(p_adjust = stats::p.adjust(p.value, method = correction_method)) |>
        dplyr::select(-contrast) |>
        tibble::as_tibble()

    # Add exposure annotations from codebook
    codebook <- MultiAssayExperiment::metadata(exposomicset)$codebook
    if (!is.null(codebook) && "variable" %in% colnames(codebook)) {
        results <- results |>
            dplyr::left_join(codebook, by = c("exposure" = "variable"))
    }

    message(sprintf(
        "Completed: %d significant associations (p_adjust < 0.05)",
        sum(results$p_adjust < 0.05, na.rm = TRUE)
    ))

    if (action == "add") {
        MultiAssayExperiment::metadata(exposomicset)$association$exposure_omics <- list(
            results_df = results,
            covariates = covariates,
            exposures = exposures,
            omics_assay = omics_assay,
            method = "limma_trend"
        )

        step_record <- list(
            run_exposure_omics_association = list(
                timestamp = Sys.time(),
                params = list(
                    exposures = exposures,
                    omics_assay = omics_assay,
                    covariates = covariates,
                    scaling_method = scaling_method,
                    correction_method = correction_method
                ),
                notes = sprintf(
                    "Tested %d exposures against %d assays using limma-trend",
                    length(exposures), length(omics_assay)
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
  covariates,
  scaling_method
) {
    exposures |>
        lapply(function(exp_var) {
            # Scale numeric exposures, keep categorical as is
            if (is.numeric(SummarizedExperiment::colData(se)[[exp_var]])) {
                SummarizedExperiment::colData(se)$.exposure <- as.numeric(
                    scale(SummarizedExperiment::colData(se)[[exp_var]])
                )
            } else {
                SummarizedExperiment::colData(se)$.exposure <-
                    SummarizedExperiment::colData(se)[[exp_var]]
            }

            formula <- if (!is.null(covariates) && length(covariates) > 0) {
                stats::as.formula(paste("~ .exposure +", paste(covariates, collapse = " + ")))
            } else {
                ~.exposure
            }

            .run_limma_trend(
                se = se,
                formula = formula,
                abundance_col = SummarizedExperiment::assayNames(se)[1],
                scaling_method = scaling_method
            ) |>
                dplyr::mutate(exposure = exp_var) |>
                dplyr::rename(
                    omics_feature = feature,
                    estimate = logFC,
                    t_statistic = t,
                    p.value = P.Value,
                    p_adjust_within = adj.P.Val
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
