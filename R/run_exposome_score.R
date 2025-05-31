#' Calculate Exposome Score for a MultiAssayExperiment
#'
#' This function calculates exposome scores using various methods (mean, median, PCA, IRT, quantile, variance)
#' based on numeric exposure variables found in the `colData` of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure data in its `colData`.
#' @param exposure_cols Optional character vector specifying which exposure columns to use. If `NULL`, all numeric columns (excluding those starting with "PC") are used.
#' @param scale Logical; if `TRUE`, standardizes (z-scores) the selected exposures before scoring. Default is `TRUE`.
#' @param score_type Character string indicating the scoring method to use. One of `"mean"`, `"median"`, `"pca"`, `"irt"`, `"quantile"`, or `"var"`.
#' @param score_column_name Optional string to rename the resulting exposome score column in the output `colData`. If `NULL`, a default name based on `score_type` is used.
#'
#' @return A `MultiAssayExperiment` object with a new column added to its `colData`, containing the computed exposome score.
#'
#' @details
#' - `"mean"`: Computes the row-wise mean across exposure values.
#' - `"median"`: Computes the row-wise median across exposure values.
#' - `"pca"`: Performs PCA on exposures and uses the first principal component.
#' - `"irt"`: Fits a unidimensional graded IRT model using decile-binned exposures.
#' - `"quantile"`: Converts each exposure to deciles (1–10) and sums them row-wise.
#' - `"var"`: Computes the row-wise variance across exposure values.
#'
#' Note: IRT-based scoring requires the `mirt` package.
#'
#' @examples
#' # Compute a PCA-based exposome score using specific exposure columns
#' # mae <- run_exposome_score(mae, exposure_cols = c("no2", "pm25"), score_type = "pca", score_column_name = "pca_score")
#'
#' @export
run_exposome_score <- function(
    expomicset,
    exposure_cols = NULL,
    scale = TRUE,
    score_type = "median",
    score_column_name = NULL){

  # Extract and preprocess colData
  message("Extracting exposure data...")
  data <- MultiAssayExperiment::colData(expomicset)  |>
    as.data.frame() |>
    dplyr::select_if(is.numeric)

  # Unselect PCs
  data <- data |>
    dplyr::select(-dplyr::starts_with("PC"))

  # Check if exposure_cols is provided
  if(!is.null(exposure_cols)){
    exposure_cols <- intersect(exposure_cols, colnames(data))

    data <- data |>
      dplyr::select(all_of(exposure_cols))
  }

  # Scale the data if user says so
  if(scale){
    data <- data |>
      dplyr::mutate_all(~ scale(.))
  }

  if(score_type == "mean"){
    # Calculate mean exposome score
    message("Calculating mean exposure scores...")

    scores <- data |>
      dplyr::mutate(exposome_score_mean = rowMeans(dplyr::across(dplyr::everything()), na.rm = TRUE)) |>
      dplyr::select(exposome_score_mean)

  } else if(score_type == "median"){
    # Calculate median exposome score
    message("Calculating median exposure scores...")

    scores <- data |>
      dplyr::mutate(exposome_score_median = matrixStats::rowMedians(
        as.matrix(
          dplyr::pick(
            dplyr::everything())),
        na.rm = TRUE))|>
      dplyr::select(exposome_score_median)

  } else if(score_type == "pca"){
    # Calculate PCA exposome score
    message("Calculating PCA exposure scores...")

    pca_result <- prcomp(data, center = TRUE)
    scores <- data.frame(pca_result$x[, 1])
    colnames(scores) <- c("exposome_score_pca")

  } else if (score_type == "irt"){
    # Calculate IRT exposome score
    message("Calculating IRT exposure scores...")

    if (!requireNamespace("mirt", quietly = TRUE)) {
      stop("Package 'mirt' is required for IRT scoring but is not installed.")
    }

    exposures_ordinal <- apply(data, 2, function(x) dplyr::ntile(x, 10))  # discretize
    model <- mirt::mirt(exposures_ordinal, 1, itemtype = "graded")
    scores <- mirt::fscores(model, full.scores.SE = FALSE)
    scores <- scores |>
      as.data.frame() |>
      setNames("exposome_score_irt")
    rownames(scores) <- rownames(data)


    } else if (score_type == "quantile"){
    # Calculate quantile exposome score
    message("Calculating quantile exposure scores...")

      quantile_data <- data |>
        dplyr::mutate(
          dplyr::across(dplyr::everything(),
                        ~ dplyr::ntile(., 10)))

      scores <- quantile_data |>
        dplyr::mutate(exposome_score_quantile = rowSums(
          dplyr::across(dplyr::everything()),
          na.rm = TRUE)) |>
        dplyr::select(exposome_score_quantile)

    }else if (score_type == "var"){
    # Calculate variance exposome score
    message("Calculating variance exposure scores...")

    scores <- data |>
      dplyr::mutate(exposome_score_var = matrixStats::rowVars(
        as.matrix(
          dplyr::pick(
            dplyr::everything())),
        na.rm = TRUE)) |>
      dplyr::select(exposome_score_var)

    } else {
    stop("Invalid score_type. Choose either 'mean', 'median', 'pca', 'irt', or 'quantile'.")
  }

  if(!is.null(score_column_name)) {
    colnames(scores) <- score_column_name
  }

  # Add scores to the MultiAssayExperiment object
  updated_col_data <- MultiAssayExperiment::colData(expomicset)  |>
    as.data.frame() |>
    # if name exists in df - overwrite the score column name
    dplyr::select(-dplyr::any_of(colnames(scores))) |>
    tibble::rownames_to_column("id") |>
    dplyr::left_join(scores |>
                       tibble::rownames_to_column("id"),
              by = "id") |>
    tibble::column_to_rownames("id")

  # Update the colData of the MultiAssayExperiment object
  MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)

  # Add analysis steps taken to metadata
  MultiAssayExperiment::metadata(expomicset)$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$steps,
    "run_exposome_score"
  )

  return(expomicset)
}


