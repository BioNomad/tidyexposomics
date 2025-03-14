#' Impute Missing Values in Exposure and Omics Data
#'
#' Imputes missing values in the `colData` (exposures) and `assays` (omics data) of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param exposure_impute_method A character string specifying the imputation method for exposure data. Options: `"median"`, `"mean"`. Default is `"median"`.
#' @param omics_impute_method A character string specifying the imputation method for omics data. Options: `"knn"`, `"mice"`, `"median"`, `"mean"`, `"dep"`, `"missforest"`. Default is `"knn"`.
#'
#' @details
#' This function:
#' - Identifies datasets with missing values using `metadata(expomicset)$na_qc`.
#' - Imputes missing exposure data using `naniar::impute_median_all()` (default) or `naniar::impute_mean_all()`.
#' - Imputes missing omics data based on `omics_impute_method`:
#'   - `"knn"`: Uses `impute::impute.knn()`.
#'   - `"mice"`: Uses `mice::mice()` with predictive mean matching.
#'   - `"median"` / `"mean"`: Uses `naniar::impute_median_all()` / `naniar::impute_mean_all()`.
#'   - `"dep"`: Uses `DEP::impute()` with `"MinProb"`, `q = 0.01` (for proteomics data).
#'   - `"missforest"`: Uses `missForest::missForest()`.
#' - Updates `colData(expomicset)` and `assays(expomicset)` with imputed values.
#'
#' @return A `MultiAssayExperiment` object with imputed missing values in exposure and omics data.
#'
#' @examples
#' \dontrun{
#' expom <- run_impute_missing(
#'   expomicset = expom,
#'   exposure_impute_method = "median",
#'   omics_impute_method = "knn"
#' )
#' }
#'
#' @export
run_impute_missing <- function(expomicset,
                           exposure_impute_method = "median",
                           omics_impute_method = "knn") {
  
  # Helper function for exposure imputation
  impute_exposure <- function(data, method) {
    numeric_data <- data |> 
      dplyr::select_if(is.numeric)
    non_numeric_data <- data |> 
      dplyr::select_if(Negate(is.numeric))
    
    if (method == "median") {
      imputed_data <- naniar::impute_median_all(numeric_data)
    } else if (method == "mean") {
      imputed_data <- naniar::impute_mean_all(numeric_data)
    } else {
      stop("Unsupported exposure impute method: ", method)
    }
    
    # Combine numeric and non-numeric data
    return(imputed_data |> 
             dplyr::bind_cols(non_numeric_data))
  }
  
  # Helper function for omics imputation
  impute_omics <- function(data, method) {
    if (method == "knn") {
      
      return(impute::impute.knn(as.matrix(data))$data |> 
               as.data.frame())
      
    } else if (method == "mice") {
      
      return(mice::mice(data,
                        m = 5, 
                        maxit = 50, 
                        method = "pmm", 
                        seed = 500))
      
    } else if (method == "median") {
      
      return(naniar::impute_median_all(data))
      
    } else if (method == "mean") {
      
      return(naniar::impute_mean_all(data))
      
    } else if (method == "dep"){
      
      return(DEP::impute(data, fun = "MinProb", q = 0.01))
      
    } else if (method == "missforest"){
      
      return(missForest::missForest(assay(data))$ximp)
      
    } else {
      stop("Unsupported omics impute method: ", method)
    }
  }
  
  # Identify datasets with missing data
  to_impute <- names(Filter(function(x) x[["all_var_sum"]] |> 
                              nrow() > 1,
                            MultiAssayExperiment::metadata(expomicset)$na_qc))
  
  omics_to_impute <- setdiff(to_impute, "exposure")
  
  # Impute exposure data if needed
  if ("exposure" %in% to_impute) {
    message("Imputing missing exposure data using ", exposure_impute_method)
    imputed_exposure <- impute_exposure(
      as.data.frame(
        MultiAssayExperiment::colData(expomicset)),
      exposure_impute_method)
    MultiAssayExperiment::colData(expomicset) <- as(imputed_exposure, "DataFrame")
  }
  
  # Impute omics data if needed
  for (omics_name in omics_to_impute) {
    message("Imputing missing omics data for ", omics_name, " using ", omics_impute_method)
    experiment <- MultiAssayExperiment::experiments(expomicset)[[omics_name]]
    assay_data <- SummarizedExperiment::assays(experiment)[[1]]
    SummarizedExperiment::assays(experiment)[[1]] <- impute_omics(
      as.data.frame(assay_data),
      omics_impute_method)
    MultiAssayExperiment::experiments(expomicset)[[omics_name]] <- experiment
  }
  
  return(expomicset)
}
