#' Run Correlation Analysis
#'
#' Computes correlations between exposures and feature types including DEGs, omics, latent factors,
#' top factor features, or principal components (PCs). Optionally computes feature–feature correlations
#' to support network analysis.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param feature_type Type of features to correlate. One of `"degs"`, `"omics"`, `"factors"`, `"factor_features"`, `"exposures"`, or `"pcs"`.
#' @param exposure_cols Optional vector of exposure column names (from `colData`) to use.
#' @param variable_map Optional mapping of features to include by assay for `omics` mode.
#' @param n_pcs Number of PCs to use when `feature_type = "pcs"`.
#' @param feature_cors Logical; if `TRUE`, compute correlations between features rather than with exposures.
#' @param robust Logical; restrict DEGs to those passing sensitivity threshold.
#' @param score_col Column name in sensitivity analysis with feature stability score.
#' @param score_thresh Threshold for filtering robust features.
#' @param correlation_method One of `"pearson"`, `"spearman"`, or `"kendall"`.
#' @param correlation_cutoff Minimum absolute correlation to retain.
#' @param cor_pval_column Column in output to filter by p-value (default: `"p.value"`).
#' @param pval_cutoff Maximum p-value or FDR threshold to retain a correlation.
#' @param deg_pval_col Column with DEG adjusted p-values.
#' @param deg_logfc_col Column with DEG log fold-changes.
#' @param deg_pval_thresh P-value cutoff for DEGs.
#' @param deg_logfc_thresh Log fold-change cutoff for DEGs.
#' @param batch_size Number of features to process per batch (default: 1500).
#' @param action Whether to `"add"` results to metadata or `"get"` as a data frame.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment` with results added to metadata.
#'         If `action = "get"`, returns a tidy `data.frame` of correlations.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_correlation(expomicset, feature_type = "omics", correlation_method = "spearman")
#' run_correlation(expomicset, feature_type = "degs", feature_cors = TRUE)
#' }
run_correlation <- function(
    expomicset,
    feature_type = c("degs", "omics", "factors", "factor_features", "exposures","pcs"),
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
    action = c("add", "get")
) {
  feature_type <- match.arg(feature_type)
  action <- match.arg(action)

  col_df <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame()
  #dplyr::select(-dplyr::starts_with("PC"))

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


  if(feature_cors &&  feature_type %in% c("degs", "omics", "factors", "factor_features")){
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
          paste0("^(", paste0(stringr::str_replace_all(exp_names, " ", "[ _]"), collapse = "|"), ")")
        ),
        # Remove matched exp_name and any underscore or space after it
        var1 = dplyr::case_when(
          !is.na(exp_name_1) ~ stringr::str_remove(var1, paste0("^", exp_name_1, "[ _]?")),
          TRUE ~ var1
        ),
        exp_name_2 = stringr::str_extract(
          var2,
          paste0("^(", paste0(stringr::str_replace_all(exp_names, " ", "[ _]"), collapse = "|"), ")")
        ),
        # Remove matched exp_name and any underscore or space after it
        var2 = dplyr::case_when(
          !is.na(exp_name_2) ~ stringr::str_remove(var2, paste0("^", exp_name_2, "[ _]?")),
          TRUE ~ var2
        )
      )

  } else{
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


  if (!feature_cors && feature_type %in% c("degs", "omics", "factors", "factor_features")) {
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
        MultiAssayExperiment::metadata(expomicset)$codebook,
        by = c("exposure"="variable")
      )

  }



  if (action == "add") {
    if (feature_cors){
      MultiAssayExperiment::metadata(expomicset)$correlation[[paste0(feature_type,"_feature_cor")]] <- correlation_df
      step_record <- list(
        run_correlation = list(
          timestamp = Sys.time(),
          params = list(
            feature_type = feature_type,
            correlation_method = correlation_method,
            correlation_cutoff = correlation_cutoff,
            pval_cutoff = pval_cutoff
          ),
          notes = paste0("Correlated ", feature_type, " features with features.")
        )
      )
      names(step_record) <- paste0("run_correlation_", feature_type)
    } else{
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

.run_correlation_omics_batches <- function(
    feature_data,
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
    cols_to_use <- intersect(features, colnames(feature_data))
    mat <- as.matrix(feature_data[, cols_to_use, drop = FALSE])
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
      dplyr::filter(abs(correlation) > correlation_cutoff) |>
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
