#' Generalized Correlation Analysis Between Exposures and Features
#'
#' Performs correlation analysis between exposures and a selected feature set (DEGs, omics, latent factors, or exposures).
#'
#' @param expomicset A MultiAssayExperiment object.
#' @param feature_type One of "degs", "omics", "factors", "factor_features", or "exposures".
#' @param exposure_cols Optional character vector of exposure variable names. If NULL, all numeric exposures are used.
#' @param variable_map Optional variable mapping table (only for omics).
#' @param robust If TRUE and feature_type == "degs", restrict to stable DEGs.
#' @param score_col Column name for stability scores in sensitivity analysis.
#' @param score_thresh Optional numeric threshold for feature stability scores.
#' @param correlation_method Correlation method ("spearman", "pearson", etc.).
#' @param correlation_cutoff Threshold for absolute correlation.
#' @param cor_pval_column Which p-value column to apply cutoff to.
#' @param pval_cutoff Significance threshold for p-values.
#' @param deg_pval_col DEG p-value column name.
#' @param deg_logfc_col DEG logFC column name.
#' @param deg_pval_thresh DEG p-value filter threshold.
#' @param deg_logfc_thresh DEG logFC filter threshold.
#' @param batch_size Features per correlation batch.
#' @param action Either "get" (return table) or "add" (update metadata).
#'
#' @return If action = "add", returns modified MultiAssayExperiment. Else, returns correlation result table.
#'
#' @export
run_correlation <- function(
    expomicset,
    feature_type = c("degs", "omics", "factors", "factor_features", "exposures"),
    exposure_cols = NULL,
    variable_map = NULL,
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
    action = c("add", "get")
) {
  feature_type <- match.arg(feature_type)
  action <- match.arg(action)

  col_df <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    dplyr::select(-dplyr::starts_with("PC"))

  exposures <- col_df |>
    dplyr::select(where(is.numeric))

  if (!is.null(exposure_cols)) {
    exposures <- exposures[, intersect(colnames(exposures), exposure_cols), drop = FALSE]
  }
  if (ncol(exposures) == 0) stop("No numeric exposures found.")

  feature_matrix <- switch(
    feature_type,
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
      variable_map),
    factors = .extract_factor_matrix(expomicset),
    exposures = .extract_exposure_matrix(
      col_df,
      exposure_cols),
    factor_features = .extract_factor_feature_matrix(expomicset)
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

  if (feature_type %in% c("degs", "omics", "factors", "factor_features")) {
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
          paste0("^(", paste0(stringr::str_replace_all(exp_names, " ", "[ _]"), collapse = "|"), ")")
        ),
        # Remove matched exp_name and any underscore or space after it
        feature = dplyr::case_when(
          !is.na(exp_name) ~ stringr::str_remove(feature, paste0("^", exp_name, "[ _]?")),
          TRUE ~ feature
        )
      )


    # Merge with exposure metadata
    correlation_df <- correlation_df |>
      dplyr::left_join(
        MultiAssayExperiment::metadata(expomicset)$var_info,
        by = c("exposure"="variable")
      )

  }



  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$correlation[[feature_type]] <- correlation_df
    step_record <- list(
      run_correlation = list(
        timestamp = Sys.time(),
        params = list(
          feature_type = feature_type,
          correlation_method = correlation_method,
          correlation_cutoff = correlation_cutoff,
          pval_cutoff = pval_cutoff
        ),
        notes = paste0("Correlated ", feature_type, " features with exposures.")
      )
    )
    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )
    return(expomicset)
  } else {
    return(correlation_df)
  }
}

.run_correlation_batches <- function(
    merged_data,
    exposure_vars,
    feature_vars,
    correlation_method,
    correlation_cutoff,
    cor_pval_column,
    pval_cutoff,
    batch_size
) {
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
                          values_to = "correlation")
    pval_df <- as.data.frame(corr_mat$P) |>
      tibble::rownames_to_column("var1") |>
      tidyr::pivot_longer(-var1,
                          names_to = "var2",
                          values_to = "p.value")

    merged <- dplyr::inner_join(corr_df,
                                pval_df,
                                by = c("var1", "var2")) |>
      dplyr::filter(var1 %in% exposure_vars,
                    var2 %in% features,
                    abs(correlation) > correlation_cutoff) |>
      dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
      dplyr::filter(!!rlang::sym(cor_pval_column) < pval_cutoff)

    correlation_results[[i]] <- merged
  }

  dplyr::bind_rows(correlation_results)
}

.extract_deg_matrix <- function(
    expomicset,
    robust,
    score_col,
    score_thresh,
    deg_pval_col,
    deg_logfc_col,
    deg_pval_thresh,
    deg_logfc_thresh) {
  da <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance
  if (is.null(da)) stop("No differential_abundance in metadata")

  da <- da |>
    dplyr::filter(
      !!rlang::sym(deg_pval_col) < deg_pval_thresh,
      abs(!!rlang::sym(deg_logfc_col)) > deg_logfc_thresh
    )

  if (robust) {
    sens <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis
    if (is.null(sens)) stop("No sensitivity_analysis in metadata")

    stable_feats <- sens$feature_stability |>
      dplyr::filter(!!dplyr::sym(score_col) > (score_thresh %||% sens$score_thresh))
    da <- dplyr::semi_join(da,
                           stable_feats,
                           by = c("exp_name", "feature"))
  }

  mats <- lapply(unique(da$exp_name), function(name) {
    se <- .update_assay_colData(expomicset, name)
    feats <- da |> dplyr::filter(exp_name == name) |> dplyr::pull(feature)
    se <- se[rownames(se) %in% feats, , drop = FALSE]
    mat <- SummarizedExperiment::assay(se) |> t() |> as.data.frame()
    colnames(mat) <- paste(name, colnames(mat), sep = "_")
    tibble::rownames_to_column(mat, "id")
  })

  purrr::reduce(mats, dplyr::full_join, by = "id")
}

.extract_omics_matrix <- function(expomicset, variable_map) {
  log2_omics <- .log2_multiassay(expomicset)
  selected <- if (!is.null(variable_map)) split(variable_map$variable, variable_map$exp_name)
  else lapply(MultiAssayExperiment::experiments(log2_omics), rownames)

  dfs <- purrr::imap(MultiAssayExperiment::experiments(log2_omics), function(se, name) {
    feats <- selected[[name]]
    if (is.null(feats)) return(NULL)
    se <- se[rownames(se) %in% feats, , drop = FALSE]
    df <- SummarizedExperiment::assay(se) |> t() |> as.data.frame()
    colnames(df) <- paste(name, colnames(df), sep = "_")
    tibble::rownames_to_column(df, "id")
  }) |> purrr::compact()

  purrr::reduce(dfs,
                dplyr::full_join,
                by = "id")
}

.extract_factor_matrix <- function(expomicset) {
  result <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results
  mat <- if (result$method == "MOFA") MOFA2::get_factors(result$result)[[1]]
  else if (result$method == "MCIA") result$result@global_scores
  else if (result$method == "MCCA") result$result@sample_scores
  else stop("Unsupported integration method")
  as.data.frame(mat) |> tibble::rownames_to_column("id")
}

.extract_factor_feature_matrix <- function(expomicset) {
  top_feats <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$top_factor_features
  if (is.null(top_feats)) stop("No top_factor_features found in metadata.")
  if (!all(c("exp_name", "feature") %in% colnames(top_feats))) {
    stop("top_factor_features must contain 'exp_name' and 'feature' columns.")
  }

  log2_omics <- .log2_multiassay(expomicset)
  selected <- split(top_feats$feature, top_feats$exp_name)

  dfs <- purrr::imap(MultiAssayExperiment::experiments(log2_omics), function(se, name) {
    feats <- selected[[name]]
    if (is.null(feats)) return(NULL)
    se <- se[rownames(se) %in% feats, , drop = FALSE]
    df <- SummarizedExperiment::assay(se) |> t() |> as.data.frame()
    colnames(df) <- paste(name, colnames(df), sep = "_")
    tibble::rownames_to_column(df, "id")
  }) |> purrr::compact()

  purrr::reduce(dfs, dplyr::full_join, by = "id")
}


.extract_exposure_matrix <- function(col_df, exposure_cols) {
  df <- col_df |> dplyr::select(where(is.numeric))
  if (!is.null(exposure_cols)) df <- df[, exposure_cols, drop = FALSE]
  tibble::rownames_to_column(df, "id")
}
