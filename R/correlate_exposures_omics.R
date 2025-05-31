#' Correlate exposures and omics features in a MultiAssayExperiment
#'
#' Performs correlation analysis between selected exposures and omics features
#' in a MultiAssayExperiment. The user can specify exactly which variables to include
#' via a mapping table to avoid name collisions and reduce computation.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param variable_map A data.frame with two columns: `variable` (variable name) and `exp_name`
#'        (either "exposure" or the name of an experiment in `experiments(expomicset)`).
#' @param correlation_method Method for correlation: "spearman" (default), "pearson", etc.
#' @param correlation_cutoff Minimum absolute correlation to keep (default 0.3).
#' @param cor_pval_column Which p-value column to use for filtering ("p.value" or "FDR").
#' @param pval_cutoff P-value threshold for significance (default 0.05).
#' @param action Whether to `"add"` the correlation results to metadata or `"get"` them (default "add").
#'
#' @return A modified MultiAssayExperiment object (if `action = "add"`), or a correlation results table (if `action = "get"`).
#'
#' @import MultiAssayExperiment SummarizedExperiment dplyr purrr tibble reshape2 Hmisc
#' @export
correlate_exposoures_omics <- function(
    expomicset,
    variable_map = NULL,
    correlation_method = "spearman",
    correlation_cutoff = 0.3,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05,
    action = "add"
) {
  message("Preparing data...")

  # Validate variable_map
  if (!is.null(variable_map)) {
    stopifnot(all(c("variable", "exp_name") %in% colnames(variable_map)))
  }

  # Extract and preprocess colData
  col_df <- MultiAssayExperiment::colData(expomicset) |> as.data.frame()

  # Select exposure variables
  if (!is.null(variable_map)) {
    exposure_vars <- variable_map$variable[variable_map$exp_name == "exposure"]
    exposure_vars <- intersect(exposure_vars, colnames(col_df))
    exposure_data <- col_df[, exposure_vars]
  } else {
    exposure_data <- dplyr::select_if(col_df, is.numeric)
    # Remove PCs
    if (any(grepl("^PC", colnames(exposure_data)))) {
      exposure_data <- dplyr::select(exposure_data, -dplyr::matches("^PC"))
    }
  }

  # Log2-transform omics
  log2_omics <- .log2_multiassay(expomicset)

  # Determine selected features per omics dataset
  if (!is.null(variable_map)) {
    omics_map <- variable_map[variable_map$exp_name != "exposure", ]
    selected_features <- split(omics_map$variable, omics_map$exp_name)
  } else {
    # Use all features if no map provided
    selected_features <- lapply(
      MultiAssayExperiment::experiments(log2_omics),
      function(se) rownames(se)
    )
  }

  # Extract and filter omics
  omics_df <- lapply(names(MultiAssayExperiment::experiments(log2_omics)), function(name) {
    se <- MultiAssayExperiment::experiments(log2_omics)[[name]]
    feats <- selected_features[[name]]
    if (is.null(feats)) return(NULL)
    se <- se[rownames(se) %in% feats, , drop = FALSE]
    df <- SummarizedExperiment::assay(se) |> t() |> as.data.frame()
    colnames(df) <- paste0(name, "_", colnames(df))  # disambiguate
    colnames(df) <- gsub(" |-", "_", colnames(df))
    df <- df |> tibble::rownames_to_column("id")
    return(df)
  }) |>
    purrr::compact() |>
    purrr::reduce(dplyr::full_join,
                  by = "id")

  # Merge with exposures
  merged_data <- exposure_data |>
    tibble::rownames_to_column("id") |>
    dplyr::left_join(omics_df, by = "id") |>
    na.omit()

  if (ncol(merged_data) > 5000) {
    message("Warning: Merged data has > 5000 variables. This may take time.")
  }

  # Correlation analysis
  message("Running correlation...")
  correlation_matrix <- merged_data |>
    dplyr::select(-id) |>
    dplyr::mutate_all(as.numeric) |>
    as.matrix() |>
    Hmisc::rcorr(type = correlation_method)

  # Melt correlation and p-values
  # correlation_df <- correlation_matrix$r |>
  #   as.data.frame() |>
  #   tibble::rownames_to_column("var1") |>
  #   reshape2::melt(id.vars = "var1") |>
  #   `colnames<-`(c("var1", "var2", "correlation"))

  correlation_df <- correlation_matrix$r |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    reshape2::melt(id.vars = "var1") |>
    `colnames<-`(c("var1", "var2", "correlation")) |>
    dplyr::mutate(
      var1 = as.character(var1),
      var2 = as.character(var2)
    )

  # pvalue_df <- correlation_matrix$P |>
  #   as.data.frame() |>
  #   tibble::rownames_to_column("var1") |>
  #   reshape2::melt(id.vars = "var1") |>
  #   `colnames<-`(c("var1", "var2", "p.value"))

  pvalue_df <- correlation_matrix$P |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    reshape2::melt(id.vars = "var1") |>
    `colnames<-`(c("var1", "var2", "p.value")) |>
    dplyr::mutate(
      var1 = as.character(var1),
      var2 = as.character(var2)
    )

  # Merge and filter results
  correlation_results <- correlation_df |>
    dplyr::inner_join(pvalue_df, by = c("var1", "var2")) |>
    dplyr::filter(abs(correlation) > correlation_cutoff) |>
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
    dplyr::filter(!!dplyr::sym(cor_pval_column) < pval_cutoff) |>
    dplyr::arrange(dplyr::desc(abs(correlation))) |>
    dplyr::mutate(
      var_min = pmin(var1, var2),
      var_max = pmax(var1, var2)
    ) |>
    dplyr::distinct(var_min, var_max, .keep_all = TRUE) |>
    dplyr::filter(var_min != var_max) |>
    # rename(var1 = var_min, var2 = var_max)
    dplyr::transmute(
      var1 = var_min,
      var2 = var_max,
      correlation,
      p.value,
      FDR
    )

  # Annotate variable types
  classify_variable <- function(varname) {
    exp_names <- names(MultiAssayExperiment::experiments(expomicset)) |>
      (\(chr)gsub(" |-", "_", chr))()

    exposures <- colnames(SummarizedExperiment::colData(expomicset))

    # Check if it's a colData variable (exposure)
    if (varname %in% exposures) {
      return("exposure")
    }

    # Check which assay prefix it came from
    for (exp_name in exp_names) {
      if (startsWith(varname, exp_name)) {
        return(exp_name)
      }
    }

    return("unknown")
  }


  correlation_results <- correlation_results |>
    dplyr::mutate(
      var1_type = purrr::map_chr(var1, classify_variable),
      var2_type = purrr::map_chr(var2, classify_variable)
    )

  correlation_results <- correlation_results |>
    dplyr::mutate(
      var1_clean = dplyr::case_when(
        var1_type != "exposure" ~ stringr::str_remove(var1, paste0("^", var1_type, "_")),
        TRUE ~ var1
      ),
      var2_clean = dplyr::case_when(
        var2_type != "exposure" ~ stringr::str_remove(var2, paste0("^", var2_type, "_")),
        TRUE ~ var2
      )
    ) |>
    dplyr::mutate(var1_type = gsub("_"," ", var1_type),
           var2_type = gsub("_"," ", var2_type)) |>
    dplyr::select(
      var1 = var1_clean,
      var2 = var2_clean,
      correlation,
      p.value,
      FDR,
      var1_type,
      var2_type
    )

  message("Correlation analysis complete.")

  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$omics_exposure_correlation <- correlation_results
    return(expomicset)
  } else if (action == "get") {
    return(correlation_results)
  } else {
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
