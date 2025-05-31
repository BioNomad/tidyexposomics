#' Extract Merged Omics and Exposure Data Frame
#'
#' This function extracts and merges exposure variables from `colData` with selected features from omics datasets
#' in a `MultiAssayExperiment` object. Optionally applies log2 transformation to omics data and restricts features based on a variable map.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param variable_map A data frame with columns `"variable"` and `"exp_name"`, indicating which variables belong to each omics or exposure domain.
#' @param log2_trans Logical; whether to log2-transform omics data. Default is `TRUE`.
#'
#' @return A data frame where rows correspond to samples, and columns contain exposure variables and log2-transformed omics features.
#'         Columns from different omics types are disambiguated using prefixes.
#'
#' @details
#' If `variable_map` is provided, it is used to select variables from both exposures and omics. If not provided, all numeric `colData` variables
#' are used as exposures (excluding variables matching `^PC`), and all omics features are included.
#'
#' @examples
#' \dontrun{
#' merged_df <- extract_omics_exposure_df(expomicset, variable_map = my_map)
#' }
#'
#' @export
extract_omics_exposure_df <- function(
    expomicset,
    variable_map,
    log2_trans=TRUE){

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

  if(log2_trans){
    # Log2-transform omics
    omics_data <- .log2_multiassay(expomicset)

  } else{
    omics_data <- expomicset
  }

  # Determine selected features per omics dataset
  if (!is.null(variable_map)) {
    omics_map <- variable_map[variable_map$exp_name != "exposure", ]
    selected_features <- split(omics_map$variable, omics_map$exp_name)
  } else {
    # Use all features if no map provided
    selected_features <- lapply(
      MultiAssayExperiment::experiments(omics_data),
      function(se) SummarizedExperiment::rownames(se)
    )
  }

  # Extract and filter omics
  omics_df <- lapply(names(MultiAssayExperiment::experiments(omics_data)), function(name) {
    se <- MultiAssayExperiment::experiments(omics_data)[[name]]
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
    na.omit() |>
    `rownames<-`(NULL) |>
    tibble::column_to_rownames("id")

  return(merged_data)

}
