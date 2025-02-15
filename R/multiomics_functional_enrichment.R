.da_multiomics_functional_enrichment <- function(
    expOmicSet, 
    pval_col = "adj.P.Val", 
    logfc_col = "logFC", 
    pval_threshold = 0.05, 
    logFC_threshold = log2(1.5), 
    ontology = "BP",  # BP (biological process), MF (molecular function), CC (cellular component)
    mirna_assays = NULL,  # User-specified miRNA assay names
    uniprot_assays = NULL  # User-specified protein assay names with UniProt IDs
) {
  require(clusterProfiler)
  require(tidyverse)
  require(org.Hs.eg.db)
  require(biomaRt)
  require(multiMiR)
  
  message("Starting multi-omics pathway enrichment analysis...")
  
  da_results <- metadata(expOmicSet)$differential_abundance
  
  if (is.null(da_results) || nrow(da_results) == 0) {
    stop("No differential abundance results found in metadata.")
  }
  
  unique_assays <- unique(da_results$assay_name)
  enrichment_results <- list()
  
  for (assay in unique_assays) {
    message("Processing assay: ", assay)
    
    # Extract all features tested for dynamic universe (excluding miRNA)
    if (is.null(mirna_assays) || !(assay %in% mirna_assays)) {
      tested_features <- rownames(experiments(expOmicSet)[[assay]])
    } else {
      tested_features <- NULL
    }
    
    # Extract significant features
    sig_features <- da_results |> 
      filter(assay_name == assay, 
             !!sym(pval_col) < pval_threshold, 
             abs(!!sym(logfc_col)) > logFC_threshold) |> 
      pull(molecular_feature) |> unique()
    
    if (length(sig_features) == 0) {
      message("No significant features for ", assay, ". Skipping...")
      next
    }
    
    # **Handle Protein Assays (Convert UniProt to Gene Symbols)**
    if (!is.null(uniprot_assays) && assay %in% uniprot_assays) {
      message("Converting UniProt IDs to gene symbols for: ", assay)
      id_map <- .convert_uniprot_to_symbol(sig_features)
      sig_features <- id_map$hgnc_symbol[!is.na(id_map$hgnc_symbol)]
      
      if (length(sig_features) == 0) {
        message("No valid gene symbols found after UniProt conversion for ", assay, ". Skipping...")
        next
      }
      # Convert full tested feature list to gene symbols
      if (!is.null(tested_features)) {
        tested_features <- .convert_uniprot_to_symbol(tested_features)$hgnc_symbol
      }
    }
    
    # **Handle miRNA Assays (Get Target Genes)**
    if (!is.null(mirna_assays) && assay %in% mirna_assays) {
      message("Retrieving target genes for miRNAs in: ", assay)
      sig_features <- .get_mirna_targets(sig_features)
      
      if (length(sig_features) == 0) {
        message("No validated miRNA targets found for ", assay, ". Skipping...")
        next
      }
    }
    
    # **Run GO Enrichment**
    message("Running GO enrichment for ", assay, " (", ontology, ")")
    enrich_res <- enrichGO(
      gene = sig_features, 
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL",
      ont = ontology, 
      universe = tested_features,  # Dynamically update universe
      pAdjustMethod = "fdr", 
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.1
    )
    
    # Store results
    enrichment_results[[assay]] <- enrich_res@result
  }
  
  # Save results
  metadata(expOmicSet)$multiomics_enrichment_results <- enrichment_results
  message("Multi-omics pathway enrichment analysis completed.")
  return(expOmicSet)
}
