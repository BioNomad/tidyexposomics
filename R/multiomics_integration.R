multiomics_integration <- function(expOmicSet,
                                   method = "MOFA",
                                   n_factors = 10,
                                   scale=TRUE) {
  require(MultiAssayExperiment)
  require(SummarizedExperiment)
  require(tidyverse)
  
  # Load method-specific packages
  if (method == "MOFA") {
    require(MOFA2)  # MOFA+ for multi-omics factor analysis
  } else if (method == "MCIA") {
    require(nipalsMCIA)  # Faster MCIA with NIPALS
  } else {
    stop("Invalid method. Choose from 'MOFA' or 'MCIA'.")
  }
  
  if(length(experiments(expOmicSet)) < 2){
    stop("Multi-Omics Integration requires at least two assays in the MultiAssayExperiment object.")
  }
  
  if(scale){
    expOmicSet <- scale_multiassay(expOmicSet)
  }
  
  message("Running multi-omics integration using ", method, "...")
  
  result <- NULL
  
  # MOFA integration
  if (method == "MOFA") {
    message("Applying MOFA+ integration...")
    
    # Create MOFA object from the MultiAssayExperiment
    mofa <- create_mofa(expOmicSet)
    
    # Set MOFA options
    model_opts <- get_default_model_options(mofa)
    model_opts$num_factors <- n_factors
    data_opts <- get_default_data_options(mofa)
    train_opts <- get_default_training_options(mofa)
    
    # Prepare & train MOFA model
    mofa <- prepare_mofa(
      object = mofa,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    outfile = file.path(tempdir(), "mofa_model.hdf5")
    mofa_trained <- run_mofa(mofa, outfile, use_basilisk = TRUE)
    
    # Load trained MOFA model
    result <- load_model(outfile)
    
    # NIPALS MCIA INTEGRATION 
  } else if (method == "MCIA") {
    message("Applying MCIA with NIPALS...")
    
    # Run NIPALS MCIA on the MultiAssayExperiment
    set.seed(42)  # Ensure reproducibility
    result <- nipals_multiblock(expOmicSet, col_preproc_method = "colprofile",
                                num_PCs = n_factors, tol = 1e-12, plots = "none")
    
  }
  
  # Store results in MultiAssayExperiment metadata
  expOmicSet@metadata$integration_results <- list(
    method = method,
    result = result
  )
  
  return(expOmicSet)
}
