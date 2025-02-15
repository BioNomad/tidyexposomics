multiomics_analysis <- function(expOmicSet, method = "MCIA", subsets = NULL, ncomp = 2, scale = TRUE) {
  library(omicade4)
  library(mogsa)
  library(mixOmics)
  library(MultiAssayExperiment)
  library(tidyverse)
  
  # Step 1: Extract omic layers from the ExperimentList
  message("Extracting omic layers from expOmicSet...")
  omics_list <- experiments(expOmicSet)
  
  if (length(omics_list) < 2) {
    stop("At least two omic layers are required for integration.")
  }
  
  # Step 2: Convert omic layers to matrices
  message("Converting omic layers to matrices...")
  omics_list <- lapply(omics_list, function(x) {
    if (inherits(x, "SummarizedExperiment")) {
      as.matrix(assay(x))
    } else if (is.data.frame(x)) {
      as.matrix(x)
    } else if (is.matrix(x)) {
      x
    } else {
      stop("Unsupported omic layer type: ", class(x))
    }
  })
  
  # Step 3: Subset omic layers (optional)
  if (!is.null(subsets)) {
    message("Sub-setting omic layers...")
    omics_list <- lapply(omics_list, function(x) {
      x[rownames(x) %in% subsets, , drop = FALSE]
    })
  }
  
  # Step 4: Validate and synchronize sample names
  message("Validating sample names across omic layers...")
  common_samples <- Reduce(intersect, lapply(omics_list, colnames))
  if (length(common_samples) == 0) {
    stop("No common samples found across omic layers.")
  }
  omics_list <- lapply(omics_list, function(x) {
    x[, common_samples, drop = FALSE]
  })
  
  # Step 5: Scale omic layers if specified
  if (scale) {
    message("Scaling omic layers...")
    omics_list <- lapply(omics_list, scale)
  }
  
  # Step 6: Perform selected method
  message("Performing multi-omic integration using ", method, "...")
  if (method == "MCIA") {
    result <- mcia(omics_list)
  } else if (method == "GCCA") {
    result <- mogsa::mbpca(omics_list, ncomp = ncomp, method = "globalScore")
  } else if (method == "PLS") {
    # if (length(omics_list) != 2) {
    #   stop("PLS requires exactly two omic layers.")
    # }
    result <- mixOmics::pls(X = omics_list, Y = colData(expOmicSet)[["pftfev1fvc_actual"]], ncomp = ncomp)
  } else {
    stop("Invalid method. Choose from 'MCIA', 'GCCA', or 'PLS'.")
  }
  
  # Step 7: Save results in metadata
  message("Saving results to expOmicSet metadata...")
  metadata(expOmicSet)$multi_omic_integration <- list(
    method = method,
    subsets = subsets,
    ncomp = ncomp,
    result = result
  )
  
  message("Multi-omic integration completed.")
  return(expOmicSet)
}
