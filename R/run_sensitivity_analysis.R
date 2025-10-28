#' Run Sensitivity Analysis for Differential Abundance
#'
#' Performs sensitivity analysis by systematically varying statistical methods,
#' scaling strategies, and model formulas (with optional bootstrap sampling) to
#' assess the stability of differentially abundant features.
#'
#' @param exposomicset A `MultiAssayExperiment` containing the assays to analyze.
#' @param base_formula The base model formula used for differential analysis.
#' @param abundance_col Character. Name of the column in the assays representing
#'  abundance. Default is `"counts"`.
#' @param methods Character vector of differential expression methods.
#' Options include `"limma_trend"` ,`"limma_voom"`, `"DESeq2"`, and `"edgeR_quasi_likelihood"`.
#' @param scaling_methods Character vector of normalization methods to try.
#'  Options include `"none"`, `"TMM"`, and `"quantile"`.
#' @param contrasts Optional list of contrasts to apply for differential testing.
#' @param covariates_to_remove Optional character vector of covariates
#' to remove from the base formula to generate model variants.
#' @param pval_col Name of the column containing p-values or adjusted
#' p-values used to define significance.
#' @param logfc_col Name of the column containing log fold changes.
#' @param pval_threshold Numeric threshold for significance. Default is 0.05.
#' @param logFC_threshold Numeric threshold for absolute log fold change.
#' Default is `log2(1)` (i.e., 0).
#' @param score_thresh Optional threshold for the selected stability metric.
#' If not provided, calculated using `score_quantile`.
#' @param score_quantile Quantile used to define the threshold
#'  if `score_thresh` is not provided. Default is 0.9.
#' @param stability_metric Character. Name of the column in
#' `feature_stability` to use as the scoring metric. Default is `"stability_score"`.
#' @param action Whether to `"add"` results to `metadata()` or
#' `"get"` them as a list. Default is `"add"`.
#' @param bootstrap_n Integer. Number of bootstrap iterations.
#' If 0, no resampling is performed. Default is 1.
#'
#' @return If `action = "add"`, returns a `MultiAssayExperiment`
#' with results stored in
#'   `metadata(exposomicset)$differential_analysis$sensitivity_analysis`.
#'   If `action = "get"`,
#'   returns a list with three elements:
#'   \describe{
#'     \item{\code{sensitivity_df}}{Data frame of all differential results
#'     across model/method combinations.}
#'     \item{\code{feature_stability}}{Data frame summarizing feature
#'     stability scores.}
#'     \item{\code{score_thresh}}{The threshold used to define stable features.}
#'   }
#'
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#'
#' # Run differential abundance
#' mae <- run_differential_abundance(
#'     exposomicset = mae,
#'     formula = ~ smoker + sex,
#'     abundance_col = "counts",
#'     method = "limma_voom",
#'     action = "add"
#' )
#'
#' # Run the sensitivity analysis
#' mae <- run_sensitivity_analysis(
#'     exposomicset = mae,
#'     base_formula = ~ smoker + sex,
#'     methods = c("limma_voom"),
#'     scaling_methods = c("none"),
#'     covariates_to_remove = "sex",
#'     pval_col = "P.Value",
#'     logfc_col = "logFC",
#'     pval_threshold = 0.05,
#'     logFC_threshold = 0,
#'     bootstrap_n = 3,
#'     action = "add"
#' )
#'
#' @export
run_sensitivity_analysis <- function(
  exposomicset,
  base_formula,
  abundance_col = "counts",
  methods = c("limma_trend", "limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
  scaling_methods = c("none", "TMM", "quantile"),
  contrasts = NULL,
  covariates_to_remove = NULL,
  pval_col = "adj.P.Val",
  logfc_col = "logFC",
  pval_threshold = 0.05,
  logFC_threshold = log2(1),
  score_thresh = NULL,
  score_quantile = 0.9,
  stability_metric = "stability_score",
  action = "add",
  bootstrap_n = 1
) {
    model_list <- .build_model_list(base_formula, covariates_to_remove)

    if (bootstrap_n > 0) {
        # Combine all raw differential results from bootstrap runs
        sensitivity_df <- purrr::map_dfr(seq_len(bootstrap_n), function(b) {
            # message(paste("Running bootstrap iteration", b, "of", bootstrap_n))
            message(sprintf("Running bootstrap iteration %d of %d", b, bootstrap_n))

            mae_b <- .resample_MAE(exposomicset)

            sensitivity_df <- .run_sensitivity_grid(
                mae_b,
                model_list,
                methods,
                scaling_methods,
                abundance_col,
                contrasts
            ) |>
                dplyr::mutate(bootstrap_id = b)
        })
    } else {
        # Run sensitivity analysis without bootstrapping
        sensitivity_df <- .run_sensitivity_grid(
            exposomicset,
            model_list,
            methods,
            scaling_methods,
            abundance_col,
            contrasts
        )
    }

    # Calculate feature stability scores
    feature_stability_df <- .calculate_feature_stability(
        sensitivity_df,
        pval_col,
        logfc_col,
        pval_threshold
    )

    # Check if stability_metric is valid
    if (!stability_metric %in% colnames(feature_stability_df)) {
        stop(sprintf("Invalid stability_metric: '%s'.", stability_metric))
    }


    # Determine score threshold if not provided
    if (is.null(score_thresh)) {
        stability_metric_col <- feature_stability_df[[stability_metric]]
        score_thresh <- quantile(
            stability_metric_col[feature_stability_df[[stability_metric]] > 0],
            score_quantile,
            na.rm = TRUE
        )
    }

    # Summarize stable features
    .summarize_stable_features(
        feature_stability_df,
        score_thresh,
        stability_metric
    )

    # Add results to MultiAssayExperiment metadata or return as list
    if (action == "add") {
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$differential_analysis$sensitivity_analysis <- list(
            sensitivity_df = sensitivity_df,
            feature_stability = feature_stability_df,
            score_thresh = score_thresh
        )
        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

        step_record <- list(
            run_sensitivity_analysis = list(
                timestamp = Sys.time(),
                params = list(
                    formula = format(base_formula),
                    methods = methods,
                    scaling_methods = scaling_methods,
                    abundance_col = abundance_col,
                    contrasts = contrasts,
                    covariates_to_remove = covariates_to_remove,
                    pval_threshold = pval_threshold,
                    logFC_threshold = logFC_threshold,
                    score_quantile = score_quantile,
                    stability_metric = stability_metric,
                    bootstrap_n = bootstrap_n
                ),
                notes = paste(
                    "Ran sensitivity analysis across",
                    bootstrap_n, "bootstrap iterations and",
                    length(methods), "methods and",
                    length(scaling_methods), "scaling strategies.",
                    if (!is.null(covariates_to_remove)) {
                        paste(
                            "Covariates removed in model variations:",
                            paste(covariates_to_remove, collapse = ", ")
                        )
                    } else {
                        "No covariates were removed from the base model."
                    }
                )
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else if (action == "get") {
        return(list(
            sensitivity_df = sensitivity_df,
            feature_stability = feature_stability_df,
            score_thresh = score_thresh
        ))
    } else {
        stop("Invalid action. Use 'add' or 'get'.")
    }
}

.run_sensitivity_grid <- function(
  exposomicset,
  model_list,
  methods,
  scalings,
  abundance_col,
  contrasts
) {
    grid <- expand.grid(
        model_name = names(model_list),
        method = methods,
        scaling = scalings,
        exp_name = names(MultiAssayExperiment::experiments(exposomicset)),
        stringsAsFactors = FALSE
    )

    results <- purrr::pmap_dfr(grid, function(model_name,
                                              method,
                                              scaling,
                                              exp_name) {
        formula <- model_list[[model_name]]
        exp <- .update_assay_colData(exposomicset, exp_name)

        if (method == "DESeq2") {
            SummarizedExperiment::assay(exp, abundance_col) <- round(
                SummarizedExperiment::assay(exp, abundance_col), 0
            )
        }

        res <- .run_da_pipeline(
            exp,
            formula,
            method,
            scaling,
            abundance_col,
            contrasts
        )

        if (is.null(res)) {
            return(NULL)
        }
        res$model <- model_name
        res$exp_name <- exp_name
        res
    })

    results <- results |>
        tibble::as_tibble()

    return(results)
}

.run_da_pipeline <- function(
  se,
  formula,
  method,
  scaling,
  abundance_col,
  contrasts
) {
    invisible(.run_se_differential_abundance(
        se = se,
        formula = formula,
        abundance_col = abundance_col,
        method = method,
        scaling_method = scaling,
        contrasts = contrasts
    ))
}

.summarize_stable_features <- function(
  feature_stability_df,
  score_thresh,
  stability_metric
) {
    sum <- feature_stability_df |>
        dplyr::group_by(exp_name) |>
        dplyr::reframe(
            n_above = sum(!!rlang::sym(stability_metric) > score_thresh),
            n = dplyr::n()
        )

    message(
        "Number Of Features Above Threshold Of ",
        round(score_thresh, 2), ":"
    )
    message("----------------------------------------")
    for (cur_exp_name in unique(feature_stability_df$exp_name)) {
        n_above <- sum |>
            dplyr::filter(exp_name == !!cur_exp_name) |>
            dplyr::pull(n_above)
        n <- sum |>
            dplyr::filter(exp_name == !!cur_exp_name) |>
            dplyr::pull(n)
        message(cur_exp_name, ": ", n_above, "/", n)
    }
}

# Resample samples in the MultiAssayExperiment object
.resample_MAE <- function(mae) {
    all_ids <- colnames(mae) |>
        unlist() |>
        unique()
    sample_ids <- sample(all_ids, replace = TRUE)
    return(mae[, sample_ids])
}

# build list of model formulas by removing covariates from base model
.build_model_list <- function(
  base_formula,
  covariates_to_remove
) {
    base_terms <- all.vars(base_formula)
    model_list <- list("full model" = base_formula)
    if (!is.null(covariates_to_remove)) {
        for (covar in covariates_to_remove) {
            reduced_terms <- setdiff(base_terms, covar)
            if (length(reduced_terms) > 1) {
                reduced_formula <- as.formula(
                    paste("~", paste(reduced_terms, collapse = " + "))
                )
                model_list[[paste("without", covar)]] <- reduced_formula
            }
        }
    }
    return(model_list)
}
