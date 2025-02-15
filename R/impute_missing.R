impute_missing <- function(expOmicSet,
                           exposure_impute_method = "median",
                           omics_impute_method = "knn") {
  require(naniar)
  require(tidyverse)
  require(impute)
  require(mice)
  
  # Helper function for exposure imputation
  impute_exposure <- function(data, method) {
    numeric_data <- select_if(data, is.numeric)
    non_numeric_data <- select_if(data, Negate(is.numeric))
    
    if (method == "median") {
      imputed_data <- impute_median_all(numeric_data)
    } else if (method == "mean") {
      imputed_data <- impute_mean_all(numeric_data)
    } else {
      stop("Unsupported exposure impute method: ", method)
    }
    
    # Combine numeric and non-numeric data
    return(bind_cols(imputed_data, non_numeric_data))
  }
  
  # Helper function for omics imputation
  impute_omics <- function(data, method) {
    if (method == "knn") {
      return(impute::impute.knn(as.matrix(data))$data %>% as.data.frame())
    } else if (method == "mice") {
      return(mice::mice(data, m = 5, maxit = 50, method = "pmm", seed = 500))
    } else if (method == "median") {
      return(impute_median_all(data))
    } else if (method == "mean") {
      return(impute_mean_all(data))
    } else if (method == "dep"){
      return(DEP::impute(data, fun = "MinProb", q = 0.01))
    } else if (method == "missforest"){
      return(missForest::missForest(assay(data))$ximp)
    } else {
      stop("Unsupported omics impute method: ", method)
    }
  }
  
  # Identify datasets with missing data
  to_impute <- names(Filter(function(x) x[["all_var_sum"]] %>% nrow() > 1, expOmicSet@metadata$na_qc))
  
  omics_to_impute <- setdiff(to_impute, "exposure")
  
  # Impute exposure data if needed
  if ("exposure" %in% to_impute) {
    message("Imputing missing exposure data using ", exposure_impute_method)
    imputed_exposure <- impute_exposure(as.data.frame(colData(expOmicSet)), exposure_impute_method)
    colData(expOmicSet) <- as(imputed_exposure, "DataFrame")
  }
  
  # Impute omics data if needed
  for (omics_name in omics_to_impute) {
    message("Imputing missing omics data for ", omics_name, " using ", omics_impute_method)
    experiment <- experiments(expOmicSet)[[omics_name]]
    assay_data <- assays(experiment)[[1]]
    assays(experiment)[[1]] <- impute_omics(as.data.frame(assay_data), omics_impute_method)
    experiments(expOmicSet)[[omics_name]] <- experiment
  }
  
  return(expOmicSet)
}
