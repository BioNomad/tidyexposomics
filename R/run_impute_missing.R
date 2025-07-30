#' Impute Missing Exposure and Omics Data in a MultiAssayExperiment
#'
#' Performs missing data imputation on both exposure variables
#' (from `colData`) and omics datasets (from `experiments`) within
#' a `MultiAssayExperiment` object.
#'
#' For exposures, numeric columns in `colData` are imputed using
#' the selected method. For omics data, assays are selected and
#' imputed individually.
#'
#' Supported imputation methods include:
#' \itemize{
#'   \item \code{"median"}: Median imputation using
#'   \code{naniar::impute_median_all}
#'   \item \code{"mean"}: Mean imputation using
#'   \code{naniar::impute_mean_all}
#'   \item \code{"knn"}: k-nearest neighbor imputation using
#'   \code{impute::impute.knn}
#'   \item \code{"mice"}: Multiple imputation using chained equations
#'   (\code{mice::mice})
#'   \item \code{"dep"}: MinProb imputation for proteomics using
#'   \code{DEP::impute}
#'   \item \code{"missforest"}: Random forest-based imputation using
#'   \code{missForest::missForest}
#'   \item \code{"lod_sqrt2"}: Substitution of missing values with
#'   LOD/sqrt(2), where LOD is the smallest non-zero value per variable
#' }
#'
#' @param expomicset A \code{MultiAssayExperiment} object containing
#' exposures and omics data.
#' @param exposure_impute_method Character. Imputation method to use for
#' exposure variables. Defaults to \code{"median"}.
#' @param exposure_cols Character vector. Names of columns in
#' \code{colData} to impute. If \code{NULL}, all numeric columns are used.
#' @param omics_impute_method Character. Imputation method to use for
#'  omics data. Defaults to \code{"knn"}.
#' @param omics_to_impute Character vector. Names of omics datasets to impute.
#' If \code{NULL}, all omics datasets are included.
#'
#' @return A \code{MultiAssayExperiment} object with imputed exposure
#' and/or omics data.
#'
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'   n_samples = 20,
#'   return_mae = TRUE
#' )
#'
#' # Introduce some missingness
#' MultiAssayExperiment::colData(mae)$exposure_pm25[sample(1:20, 5)] <- NA
#'
#' # Filter features and exposures with high missingness
#' mae <- run_impute_missing(
#'   expomicset = mae,
#'   exposure_impute_method = "median"
#' )
#'
#' @export
run_impute_missing <- function(expomicset,
                               exposure_impute_method = "median",
                               exposure_cols = NULL,
                               omics_impute_method = NULL,
                               omics_to_impute = NULL) {

  # Helper for LOD/sqrt(2) imputation
  impute_lod_sqrt2 <- function(data) {
    data[] <- lapply(data, function(col) {
      if (is.numeric(col)) {
        lod <- min(col[col > 0], na.rm = TRUE)
        imputed <- ifelse(is.na(col), lod / sqrt(2), col)
        return(imputed)
      } else {
        return(col)
      }
    })
    return(as.data.frame(data))
  }

  # General imputation helper
  impute_data <- function(data, method) {
    if (method == "median") {
      return(naniar::impute_median_all(data))
    } else if (method == "mean") {
      return(naniar::impute_mean_all(data))
    } else if (method == "knn") {
      return(as.data.frame(impute::impute.knn(as.matrix(data))$data))
    } else if (method == "mice") {
      return(mice::complete(mice::mice(data, m = 5, maxit = 50,
                                       method = "pmm", seed = 500)))
    } else if (method == "dep") {
      return(DEP::impute(data, fun = "MinProb", q = 0.01))
    } else if (method == "missforest") {
      return(missForest::missForest(as.matrix(data))$ximp)
    } else if (method == "lod_sqrt2") {
      return(impute_lod_sqrt2(data))
    } else {
      stop("Unsupported imputation method: ", method)
    }
  }

  # Impute selected exposure columns
  if ("exposure" %in% names(
    MultiAssayExperiment::metadata(expomicset)$quality_control$na_qc)) {
    message("Imputing exposure data using method: ", exposure_impute_method)

    metadata_df <- as.data.frame(MultiAssayExperiment::colData(expomicset))
    # if (is.null(exposure_cols)) {
    #   exposure_cols <- names(metadata_df)[sapply(metadata_df, is.numeric)]
    # }
    if (is.null(exposure_cols)) {
      exposure_cols <- names(metadata_df)[
        vapply(metadata_df, is.numeric, FUN.VALUE = logical(1))
      ]
    }


    data_to_impute <- metadata_df[, exposure_cols, drop = FALSE]
    imputed <- impute_data(data_to_impute, exposure_impute_method)

    metadata_df[, exposure_cols] <- imputed
    MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(
      metadata_df)
  }

  if(!is.null(omics_to_impute)){
    # Impute omics data
    all_omics <- setdiff(names(
      MultiAssayExperiment::experiments(expomicset)), "exposure")
    omics_to_use <- if (is.null(omics_to_impute)){
      all_omics
      } else {
        intersect(all_omics, omics_to_impute)
      }

    for (omics_name in omics_to_use) {
      message("Imputing omics dataset: ",
              omics_name,
              " using method: ",
              omics_impute_method)
      experiment <- MultiAssayExperiment::experiments(expomicset)[[omics_name]]
      assay_data <- as.data.frame(SummarizedExperiment::assays(experiment)[[1]])
      imputed_data <- impute_data(assay_data, omics_impute_method)
      SummarizedExperiment::assays(experiment)[[1]] <- as.matrix(imputed_data)
      MultiAssayExperiment::experiments(expomicset)[[omics_name]] <- experiment
    }
  }

  # Add analysis steps taken to metadata
  step_record <- list(run_impute_missing=list(
    timestamp = Sys.time(),
    params = list(exposure_impute_method = exposure_impute_method,
                  exposure_cols = exposure_cols,
                  omics_impute_method = omics_impute_method,
                  omics_to_impute = omics_to_impute),
    notes = ""))

  MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$summary$steps,
    step_record
  )

  return(expomicset)
}

