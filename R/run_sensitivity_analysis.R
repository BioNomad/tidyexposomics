#' Perform Sensitivity Analysis for Differential Abundance
#'
#' Runs differential abundance testing across multiple models, statistical methods,
#' scaling approaches, and filtering criteria to assess feature stability.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param base_formula A formula specifying the base model for differential abundance analysis.
#' @param abundance_col A character string specifying the assay column to use for abundance values. Default is `"counts"`.
#' @param methods A character vector of methods for differential abundance testing (e.g., `"limma_voom"`, `"DESeq2"`, `"edgeR_quasi_likelihood"`).
#' @param scaling_methods A character vector of normalization methods to apply (e.g., `"none"`, `"TMM"`, `"quantile"`).
#' @param min_counts_range A numeric vector of minimum count thresholds to test. Default is `c(1, 5, 10)`.
#' @param min_proportion_range A numeric vector of minimum sample proportion thresholds to test. Default is `c(0.1, 0.2, 0.3)`.
#' @param contrasts A character vector specifying contrasts for the differential analysis. Default is `NULL`.
#' @param covariates_to_remove A character vector of covariates to iteratively remove from the model for sensitivity testing. Default is `NULL`.
#' @param pval_col Name of the column containing adjusted p-values. Default is `"adj.P.Val"`.
#' @param logfc_col Name of the column containing log fold changes. Default is `"logFC"`.
#' @param pval_threshold Numeric threshold to consider a p-value significant. Default is `0.05`.
#' @param logFC_threshold Minimum absolute log fold change for filtering. Default is `log2(1.5)`.
#' @param score_thresh A numeric value specifying a fixed threshold for feature stability. If `NULL`, it is determined based on `score_quantile`.
#' @param score_quantile A numeric value (between 0 and 1) specifying the quantile to define stability threshold. Default is `0.9`.
#' @param stability_metric The name of the column in the stability table to use for scoring (e.g., `"stability_score"`, `"logp_weighted_score"`). Default is `"stability_score"`.
#' @param action Whether to `"add"` results to the `MultiAssayExperiment` metadata or `"get"` them as a list. Default is `"add"`.
#' @param bootstrap_n Number of bootstrap iterations. Set to `0` to disable bootstrapping. Default is `0`.
#'
#' @return Either a `MultiAssayExperiment` object with results added to metadata (if `action = "add"`) or a list with:
#' \item{sensitivity_df}{A dataframe of differential abundance results across all tested conditions.}
#' \item{feature_stability}{A dataframe summarizing stability scores for each feature.}
#' \item{score_thresh}{The threshold used for stable feature selection.}
#'
#' @details
#' The stability score can be calculated in various ways. By default, `stability_score` is used, which multiplies the proportion of significant p-values across runs (`presence_rate`) by a consistency metric (`effect_consistency`). Users may specify alternative metrics like `logp_weighted_score` for more conservative feature prioritization.
#'
#' If `bootstrap_n > 0`, the entire sensitivity analysis is repeated `bootstrap_n` times on resampled versions of the dataset, and summary statistics are returned for each feature's stability.
#'
#' @examples
#' \dontrun{
#' result <- run_sensitivity_analysis(
#'   expomicset = expom,
#'   base_formula = ~ condition + batch,
#'   methods = c("limma_voom", "DESeq2"),
#'   min_counts_range = c(5, 10),
#'   stability_metric = "logp_weighted_score",
#'   action = "get"
#' )
#' }
#'
#' @export

run_sensitivity_analysis <- function(
    expomicset,
    base_formula,
    abundance_col = "counts",
    methods = c("limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
    scaling_methods = c("none", "TMM", "quantile"),
    min_counts_range = c(1, 5, 10),
    min_proportion_range = c(0.1, 0.2, 0.3),
    contrasts = NULL,
    covariates_to_remove = NULL,
    pval_col = "adj.P.Val",
    logfc_col = "logFC",
    pval_threshold = 0.05,
    logFC_threshold = log2(1.5),
    score_thresh = NULL,
    score_quantile = 0.9,
    stability_metric = "stability_score",
    action = "add",
    bootstrap_n = 0
) {
  if (bootstrap_n > 0) {
    return(.bootstrap_analysis(
      expomicset = expomicset,
      base_formula = base_formula,
      abundance_col = abundance_col,
      methods = methods,
      scaling_methods = scaling_methods,
      min_counts_range = min_counts_range,
      min_proportion_range = min_proportion_range,
      contrasts = contrasts,
      covariates_to_remove = covariates_to_remove,
      pval_col = pval_col,
      logfc_col = logfc_col,
      pval_threshold = pval_threshold,
      logFC_threshold = logFC_threshold,
      score_thresh = score_thresh,
      score_quantile = score_quantile,
      bootstrap_n = bootstrap_n,
      action = action
    ))
  }

  model_list <- .build_model_list(base_formula,
                                  covariates_to_remove)
  sensitivity_df <- .run_sensitivity_grid(
    expomicset,
    model_list,
    methods,
    scaling_methods,
    min_counts_range,
    min_proportion_range,
    abundance_col,
    contrasts
  )

  feature_stability_df <- .calculate_feature_stability(
    sensitivity_df,
    pval_col,
    logfc_col,
    pval_threshold
  )

  if (!stability_metric %in% colnames(feature_stability_df)) {
    stop(paste0("Invalid stability_metric: '", stability_metric, "'. Must be one of: ", paste(names(feature_stability_df), collapse = ", ")))
  }

  if (is.null(score_thresh)) {
    score_thresh <- quantile(feature_stability_df[[stability_metric]], score_quantile)
  }

  .summarize_stable_features(feature_stability_df, score_thresh, stability_metric)

  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis <- list(
      sensitivity_df = sensitivity_df,
      feature_stability = feature_stability_df,
      score_thresh = score_thresh
    )

    # Add step record
    step_record <- list(
      run_sensitivity_analysis = list(
        timestamp = Sys.time(),
        params = list(
          formula = format(base_formula),
          methods = methods,
          scaling_methods = scaling_methods,
          abundance_col = abundance_col,
          min_counts_range = min_counts_range,
          min_proportion_range = min_proportion_range,
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
          length(methods), "methods and",
          length(scaling_methods), "scaling strategies, testing",
          length(min_counts_range), "minimum count thresholds and",
          length(min_proportion_range), "sample proportion thresholds.",
          if (!is.null(covariates_to_remove)) {
            paste("Covariates removed in model variations:", paste(covariates_to_remove, collapse = ", "))
          } else {
            "No covariates were removed from the base model."
          }
        )
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)
  }
  else if (action == "get") {
    return(list(
      sensitivity_df = sensitivity_df,
      feature_stability = feature_stability_df,
      score_thresh = score_thresh
    ))
  } else {
    stop("invalid action. use 'add' or 'get'.")
  }
}

# build list of model formulas by removing covariates from base model
.build_model_list <- function(base_formula,
                              covariates_to_remove) {
  base_terms <- all.vars(base_formula)
  model_list <- list("full model" = base_formula)
  if (!is.null(covariates_to_remove)) {
    for (covar in covariates_to_remove) {
      reduced_terms <- setdiff(base_terms, covar)
      if (length(reduced_terms) > 1) {
        reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
        model_list[[paste("without", covar)]] <- reduced_formula
      }
    }
  }
  return(model_list)
}

# get features passing abundance thresholds
.get_abundant_features <- function(se, min_counts, min_prop) {
  tidybulk::identify_abundant(se,
                              minimum_counts = min_counts,
                              minimum_proportion = min_prop)
}

# run differential abundance pipeline for one configuration
.run_da_pipeline <- function(
    se,
    formula,
    method,
    scaling,
    min_counts,
    min_prop,
    abundance_col,
    contrasts) {
  .run_se_differential_abundance(
    se = se,
    formula = formula,
    abundance_col = abundance_col,
    method = method,
    scaling_method = scaling,
    min_counts = min_counts,
    min_proportion = min_prop,
    contrasts = contrasts
  )
}

# run sensitivity grid across all configurations using pmap
.run_sensitivity_grid <- function(
    expomicset,
    model_list,
    methods,
    scalings,
    min_counts_range,
    min_prop_range,
    abundance_col,
    contrasts
) {
  grid <- expand.grid(
    model_name = names(model_list),
    method = methods,
    scaling = scalings,
    min_counts = min_counts_range,
    min_prop = min_prop_range,
    exp_name = names(MultiAssayExperiment::experiments(expomicset)),
    stringsAsFactors = FALSE
  )

  results <- purrr::pmap_dfr(grid, function(
    model_name,
    method,
    scaling,
    min_counts,
    min_prop,
    exp_name) {
    formula <- model_list[[model_name]]
    exp <- .update_assay_colData(expomicset, exp_name)
    abundant <- .get_abundant_features(exp, min_counts, min_prop)
    n_features <- sum(S4Vectors::elementMetadata(abundant)[[".abundant"]], na.rm = TRUE)
    if (n_features < 2) return(NULL)
    if (method == "DESeq2") {
      SummarizedExperiment::assay(exp, abundance_col) <- round(SummarizedExperiment::assay(exp, abundance_col), 0)
    }
    res <- .run_da_pipeline(
      exp,
      formula,
      method,
      scaling,
      min_counts,
      min_prop,
      abundance_col,
      contrasts)
    if (is.null(res)) return(NULL)
    res$model <- model_name
    res$exp_name <- exp_name
    res
  })

  return(results)
}

# summarize and report stable features above threshold
.summarize_stable_features <- function(
    feature_stability_df,
    score_thresh,
    stability_metric) {
  sum <- feature_stability_df |>
    dplyr::group_by(exp_name) |>
    dplyr::reframe(
      n_above = sum(!!rlang::sym(stability_metric) > score_thresh),
      n = dplyr::n()
    )

  message("Number Of Features Above Threshold Of ", round(score_thresh, 2), ":")
  message("----------------------------------------")
  for (exp_name in unique(feature_stability_df$exp_name)) {
    n_above <- sum |>
      dplyr::filter(exp_name == !!exp_name) |>
      dplyr::pull(n_above)
    n <- sum |>
      dplyr::filter(exp_name == !!exp_name) |>
      dplyr::pull(n)
    message(exp_name, ": ", n_above, "/", n)
  }
}

# run sensitivity analysis with bootstrapping
.bootstrap_analysis <- function(
    expomicset,
    base_formula,
    abundance_col,
    methods,
    scaling_methods,
    min_counts_range,
    min_proportion_range,
    contrasts,
    covariates_to_remove,
    pval_col,
    logfc_col,
    pval_threshold,
    logFC_threshold,
    score_thresh,
    score_quantile,
    bootstrap_n,
    action
) {
  message("Running Bootstrapped Sensitivity Analysis With ", bootstrap_n, " iterations...")

  boot_results <- purrr::map_dfr(1:bootstrap_n, function(b) {
    mae_b <- .resample_MAE(expomicset)
    res_b <- run_sensitivity_analysis(
      expomicset = mae_b,
      base_formula = base_formula,
      abundance_col = abundance_col,
      methods = methods,
      scaling_methods = scaling_methods,
      min_counts_range = min_counts_range,
      min_proportion_range = min_proportion_range,
      contrasts = contrasts,
      covariates_to_remove = covariates_to_remove,
      pval_col = pval_col,
      logfc_col = logfc_col,
      pval_threshold = pval_threshold,
      logFC_threshold = logFC_threshold,
      score_thresh = score_thresh,
      score_quantile = score_quantile,
      action = "get"
    )
    res_b$feature_stability |>
      dplyr::mutate(bootstrap_id = b)
  })

  # boot_summary <- boot_results |>
  #   dplyr::group_by(feature, exp_name) |>
  #   dplyr::summarise(
  #     mean_stability = mean(stability_score, na.rm = TRUE),
  #     freq_selected = mean(stability_score > score_thresh, na.rm = TRUE),
  #     sd_stability = sd(stability_score, na.rm = TRUE),
  #     .groups = "drop"
  #   )

  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis <- list(
      feature_stability = boot_results,
      score_thresh = score_thresh
    )

    # Add step record
    step_record <- list(
      run_sensitivity_analysis_bootstrap = list(
        timestamp = Sys.time(),
        params = list(
          formula = format(base_formula),
          methods = methods,
          scaling_methods = scaling_methods,
          abundance_col = abundance_col,
          min_counts_range = min_counts_range,
          min_proportion_range = min_proportion_range,
          contrasts = contrasts,
          covariates_to_remove = covariates_to_remove,
          pval_threshold = pval_threshold,
          logFC_threshold = logFC_threshold,
          score_quantile = score_quantile,
          stability_metric = "stability_score",
          bootstrap_n = bootstrap_n
        ),
        notes = paste(
          "Performed bootstrapped sensitivity analysis (", bootstrap_n, " iterations) using",
          length(methods), "methods and",
          length(scaling_methods), "scaling strategies, with feature-level stability scores recorded per run.",
          if (!is.null(covariates_to_remove)) {
            paste("Covariates removed in model variations:", paste(covariates_to_remove, collapse = ", "))
          } else {
            "No covariates were removed from the base model."
          }
        )
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)
  }
  else if (action == "get") {
    return(list(
      feature_stability = boot_results,
      score_thresh = score_thresh
    ))
  } else {
    stop("invalid action. use 'add' or 'get'.")
  }
}


# resample columns (samples) of multiassayexperiment object
.resample_MAE <- function(mae) {
  sample_ids <- sample(colnames(mae), replace = TRUE)
  return(mae[, sample_ids])
}
