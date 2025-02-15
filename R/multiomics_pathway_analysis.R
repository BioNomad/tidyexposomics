multiomics_pathway_analysis <- function(expOmicSet,
                                        omic, 
                                        factors = 1:3, 
                                        path.database = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/c2.cp.reactome.v6.2.symbols.gmt") {
  require(fgsea)
  require(stringr)
  
  message("Extracting feature loadings for pathway analysis...")
  
  # Load GSEA pathway database
  pathways <- fgsea::gmtPathways(gmt.file = path.database)
  
  # Get the integration results
  integration_results <- expOmicSet@metadata$integration_results
  method <- integration_results$method
  
  # Extract feature loadings based on integration method
  if (method == "MOFA") {
    message("Detected MOFA+ results, extracting weights...")
    feature_loadings <- get_weights(integration_results$result, view = omic)
  } else if (method == "MCIA") {
    message("Detected MCIA results, extracting block loadings...")
    feature_loadings <- integration_results$result@block_loadings[[omic]]
  } else {
    stop("Unsupported integration method: ", method)
  }
  
  if (is.null(feature_loadings)) {
    stop("No loadings found for the specified omic: ", omic)
  }
  
  # Clean feature names (if necessary)
  rownames(feature_loadings) <- str_remove(rownames(feature_loadings), "_[0-9]*_.*")
  
  # Run GSEA on each factor
  gsea_results <- lapply(factors, function(factor) {
    fgsea::fgseaMultilevel(
      pathways = pathways, 
      stats = feature_loadings[, factor],  # Use loadings as ranking scores
      nPermSimple = 10000,
      minSize = 15,
      nproc = 2, 
      maxSize = 500
    )
  })
  
  names(gsea_results) <- paste0("Factor_", factors)
  
  # Store results in metadata
  expOmicSet@metadata$pathway_analysis <- gsea_results
  
  return(expOmicSet)
}
