scale_multiassay <- function(expOmicSet) {
  require(MultiAssayExperiment)
  require(SummarizedExperiment)
  require(tidyverse)
  
  message("Scaling each assay in MultiAssayExperiment...")
  
  # Apply scaling to each assay
  scaled_experiments <- lapply(experiments(expOmicSet), function(assay_obj) {
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
  scaled_expOmicSet <- MultiAssayExperiment(
    experiments = scaled_experiments, 
    colData = colData(expOmicSet),
    metadata = metadata(expOmicSet))
  
  return(scaled_expOmicSet)
}