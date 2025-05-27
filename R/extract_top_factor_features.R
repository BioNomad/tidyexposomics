#' Extract Top Contributing Features for Factors
#'
#' Identifies the most influential features for specified factors using MOFA+ or MCIA
#' integration results. Features are selected based on either a percentile cutoff
#' or an absolute loading threshold.
#'
#' @param expomicset A `MultiAssayExperiment` object containing integration results.
#' @param factors A character vector specifying the factors of interest.
#' @param method A character string specifying the feature selection method (`"percentile"` or `"threshold"`). Default is `"percentile"`.
#' @param percentile A numeric value between 0 and 1 indicating the percentile threshold for feature selection when `method = "percentile"`. Default is `0.9`.
#' @param threshold A numeric value specifying the absolute loading cutoff for feature selection when `method = "threshold"`. Default is `0.3`.
#' @param action A character string indicating whether to return results (`"get"`) or add them to metadata (`"add"`). Default is `"add"`.
#'
#' @details
#' The function extracts factor loadings from `metadata(expomicset)`, applies filtering based on
#' the selected method, and identifies top contributing features for each specified factor.
#' Features can be selected using:
#' - **Percentile-based filtering** (`method = "percentile"`): Selects features with absolute loadings above a specified percentile.
#' - **Threshold-based filtering** (`method = "threshold"`): Selects features with absolute loadings exceeding a fixed value.
#'
#' @return If `action = "add"`, returns the modified `expomicset` with selected features stored in metadata.
#' If `action = "get"`, returns a data frame containing:
#' \item{feature}{The selected feature contributing to the factor.}
#' \item{factor}{The factor to which the feature contributes.}
#' \item{loading}{The factor loading value of the feature.}
#' \item{exp_name}{The experiment from which the feature originated.}
#'
#' @examples
#' \dontrun{
#' results <- extract_top_factor_features(
#'   expomicset = expom,
#'   factors = c("Factor1", "Factor2"),
#'   method = "percentile",
#'   percentile = 0.9,
#'   action = "get"
#' )
#' }
#'
#' @export
extract_top_factor_features <- function(
    expomicset,
    factors,
    method = "percentile",
    percentile = 0.9,
    threshold = 0.3,
    action = "add") {

  message("Extracting top contributing features for specified factors...")

  # Get integration results
  integration_results <- MultiAssayExperiment::metadata(expomicset)$integration_results
  method_used <- integration_results$method

  # Extract factor loadings based on method
  if (method_used == "MOFA") {
    message("Using MOFA+ factor loadings...")
    loadings <- MOFA2::get_weights(integration_results$result)
  } else if (method_used == "MCIA") {
    message("Using MCIA block loadings...")
    loadings <- integration_results$result@block_loadings
  } else {
    stop("Unsupported integration method: ", method_used)
  }

  # Convert factor loadings to long format
  loadings_df <- loadings |>
    purrr::map2(names(loadings),
                \(df, exp_name) {
      df <- as.data.frame(df)
      df$exp_name <- exp_name
      df$feature <- rownames(df)
      rownames(df) <- NULL
      return(df)
    }) |>
    dplyr::bind_rows() |>
    tidyr::pivot_longer(
      -c(feature, exp_name),
      names_to = "factor",
      values_to = "loading")

  # Ensure factor names are treated as characters and match MCIA format
  factors <- as.character(factors)
  loadings_df <- loadings_df |>
    dplyr::mutate(factor = as.character(factor))

  # Filter for user-specified factors
  loadings_df <- loadings_df |>
    dplyr::filter(factor %in% factors)

  # Apply thresholding method
  # if (method == "percentile") {
  #   message(
  #     "Applying percentile-based filtering (",
  #     percentile * 100,
  #     "%)...")
  #
  #   factor_thresholds <- loadings_df |>
  #     dplyr::group_by(factor) |>
  #     dplyr::summarize(
  #       threshold = quantile(abs(loading),
  #                            percentile,
  #                            na.rm = TRUE),
  #               .groups = "drop")
  # } else if (method == "threshold") {
  #   message(
  #     "Applying raw threshold-based filtering (>|",
  #     threshold,
  #     "|)...")
  #   factor_thresholds <- tibble::tibble(
  #     factor = unique(loadings_df$factor),
  #     threshold = threshold)
  # } else {
  #   stop("Invalid method. Choose 'percentile' or 'threshold'.")
  # }

  # Merge computed thresholds with loadings
  # loadings_df <- loadings_df |>
  #   dplyr::left_join(
  #   factor_thresholds,
  #   by = "factor")

  # Apply filtering
  filtered_features <- loadings_df |>
    #dplyr::filter(abs(loading) > threshold) |>
    dplyr::select(feature,
                  factor,
                  loading,
                  exp_name) |>
    dplyr::distinct() |>
    group_by(factor) |>
    mutate(rank=percent_rank(abs(loading))) |>
    arrange(desc(rank)) |>
    filter(rank>percentile)

  message("Selected ",
          nrow(filtered_features),
          " features contributing to specified factors.")

  if(action=="add"){
    # Store selected features
    MultiAssayExperiment::metadata(expomicset)$top_factor_features <- filtered_features
    return(expomicset)
  }else if (action=="get"){
    return(filtered_features)
  }else{
    stop("Invalid action. Choose 'add' or 'get'.")
  }
}
