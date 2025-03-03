scale_multiassay <- function(expomicset) {
  require(MultiAssayExperiment)
  require(SummarizedExperiment)
  require(tidyverse)
  
  message("Scaling each assay in MultiAssayExperiment...")
  
  # Apply scaling to each assay
  scaled_experiments <- lapply(experiments(expomicset), function(assay_obj) {
    if (inherits(assay_obj, "SummarizedExperiment")) {
      assay_mat <- assay(assay_obj)
      scaled_mat <- scale(assay_mat)  # Standardize (Z-score)
      assay(assay_obj) <- scaled_mat
      return(assay_obj)
    } else if (is.matrix(assay_obj)) {
      return(scale(assay_obj))  # Directly scale matrices
    } else {
      stop("Unsupported assay type. Only SummarizedExperiment and matrices are supported.")
    }
  })
  
  # Create a new MultiAssayExperiment with scaled data
  scaled_expomicset <- MultiAssayExperiment(
    experiments = scaled_experiments, 
    colData = colData(expomicset),
    metadata = metadata(expomicset))
  
  return(scaled_expomicset)
}