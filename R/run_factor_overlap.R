#' Identify and Annotate Shared Top Features Across Integration Factors
#'
#' Identifies top features shared across factors based on integration method.
#' For MOFA/MCIA, takes intersection across factors. For DIABLO/RGCCA,
#' takes features recurring in ≥2 block-specific components.
#'
#' @param expomicset A `MultiAssayExperiment` with integration results
#' and top factor features.
#' @param robust Logical; if `TRUE`, uses sensitivity score. Otherwise,
#'  uses DEG thresholds.
#' @param stability_score Optional numeric threshold
#' (overrides default from metadata).
#' @param score_col Column name for sensitivity score.
#' Default is `"stability_score"`.
#' @param pval_thresh DEG p-value threshold (if `robust = FALSE`).
#' Default is `0.05`.
#' @param logfc_thresh DEG logFC threshold (if `robust = FALSE`).
#'  Default is `log2(1.5)`.
#' @param pval_col Column name for p-value. Default is `"padj"`.
#' @param logfc_col Column name for logFC. Default is `"logFC"`.
#' @param action `"add"` to return modified object, `"get"` to
#' return data.frame.
#'
#' @return Modified `MultiAssayExperiment` or `data.frame` of
#' shared top features.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'    n_samples = 20,
#'    return_mae=TRUE
#'   )
#'
#' # perform multiomics integration
#' mae <- run_multiomics_integration(
#'       mae,
#'       method = "MCIA",
#'       n_factors = 3)
#'
#'
#' # identify the features that contribute most to the factors
#' mae <- extract_top_factor_features(
#'   mae,
#'   factors = c("V1","V2","V3"),
#'   method = "percentile",
#'   percentile = 0.5,
#'   action = "add"
#' )
#'
#' # perform differential abundance analysis
#' mae <- run_differential_abundance(
#'   expomicset = mae,
#'   formula = ~ smoker + sex,
#'   abundance_col = "counts",
#'   method = "limma_voom",
#'   action = "add"
#' )
#'
#' # determine the overlap in features
#' mae <- mae |>
#'   run_factor_overlap(
#'   robust = FALSE,
#'   pval_col = "adj.P.Val")
#'
#' @export
run_factor_overlap <- function(
    expomicset,
    robust = TRUE,
    stability_score = NULL,
    score_col = "stability_score",
    pval_thresh = 0.05,
    logfc_thresh = log2(1.5),
    pval_col = "padj",
    logfc_col = "logFC",
    action = "add") {

  # Ensure top features exist
  if (!"top_factor_features" %in% names(MultiAssayExperiment::metadata(
    expomicset)$multiomics_integration)) {
    stop("Please run 'extract_top_factor_features()' first.")
  }

  method <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method

  top_factor_features <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$top_factor_features |>
    dplyr::mutate(exp_name_feature = paste(exp_name, feature, sep = "_"))

  # Identify common features across factors
  common_features <- if (method %in% c("DIABLO", "RGCCA")) {
    top_factor_features |>
      dplyr::group_by(exp_name, feature) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n > 1) |>
      dplyr::mutate(exp_name_feature = paste(exp_name, feature, sep = "_")) |>
      dplyr::pull(exp_name_feature)
  } else {
    top_factor_features |>
      (\(df) split(df, df$factor))() |>
      purrr::map(~dplyr::pull(.x, exp_name_feature)) |>
      purrr::reduce(intersect)
  }

  top_factor_features <- top_factor_features |>
    dplyr::filter(exp_name_feature %in% common_features)

  # Annotate with DEG/robust marker
  da_res <- if (robust) {
    score <- stability_score %||%
      purrr::pluck(MultiAssayExperiment::metadata(expomicset),
                   "differential_analysis",
                   "sensitivity_analysis",
                   "score_thresh")

    purrr::pluck(MultiAssayExperiment::metadata(expomicset),
                 "differential_analysis",
                 "sensitivity_analysis",
                 "feature_stability") |>
      dplyr::filter(!!rlang::sym(score_col) > score) |>
      dplyr::mutate(exp_name_feature = paste(exp_name, feature, sep = "_"))
  } else {
    purrr::pluck(MultiAssayExperiment::metadata(expomicset),
                 "differential_analysis",
                 "differential_abundance") |>
      dplyr::filter(
        !!rlang::sym(pval_col) < pval_thresh,
        abs(!!rlang::sym(logfc_col)) > logfc_thresh
      ) |>
      dplyr::mutate(exp_name_feature = paste(exp_name, feature, sep = "_"))
  }

  top_factor_features <- top_factor_features |>
    dplyr::mutate(is_deg = exp_name_feature %in% da_res$exp_name_feature)

  # Optional: Add annotations using pivot_feature()
  top_factor_features <- top_factor_features |>
    dplyr::left_join(
      pivot_feature(expomicset),
      by = c("feature" = ".feature",
             "exp_name" = ".exp_name")
    )

  message(glue::glue(
    "Found {length(unique(top_factor_features$exp_name_feature))}
    common features across factors."))

  if (action == "add") {
    all_metadata <- MultiAssayExperiment::metadata(expomicset)
    all_metadata$multiomics_integration$common_top_factor_features <-
      top_factor_features
    MultiAssayExperiment::metadata(expomicset) <- all_metadata

    step_record <- list(
      run_factor_overlap = list(
        timestamp = Sys.time(),
        params = list(
          robust = robust,
          stability_score = stability_score,
          score_col = score_col,
          pval_thresh = pval_thresh,
          logfc_thresh = logfc_thresh
        ),
        notes = if (robust) {
          paste0("Annotated based on '",
                 score_col,
                 "' score > ",
                 round(score, 3), ".")
        } else {
          paste0("Annotated using p<",
                 pval_thresh,
                 " and |logFC|>",
                 round(logfc_thresh, 3), ".")
        }
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)
  } else {
    return(top_factor_features)
  }
}


