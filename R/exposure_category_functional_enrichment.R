exposure_category_functional_enrichment <- function(
    expOmicSet,
    cor_df = "degs", # Can be "degs" or "factors"
    pval_col = "p.value",
    cor_threshold = 0.3,  # Absolute correlation cutoff
    pval_threshold = 0.05,  
    ontology = "BP",  # BP (biological process), MF (molecular function), CC (cellular component)
    mirna_assays = NULL,  # User-specified miRNA assay names
    uniprot_assays = NULL,  # User-specified protein assay names with UniProt IDs
    universe_background = TRUE,  # Whether to use a background gene universe (disabled for miRNA)
    orgdb = org.Hs.eg.db,  # Organism annotation database
    keytype = "SYMBOL",  # Key type for enrichment analysis
    p_adjust_method = "fdr",  # Method for multiple testing correction
    pvalue_cutoff = 0.05,  # p-value threshold for enrichment
    qvalue_cutoff = 0.1  # q-value threshold for enrichment
) {
  require(clusterProfiler)
  require(tidyverse)
  require(org.Hs.eg.db)
  require(biomaRt)
  require(multiMiR)
  
  message("Starting enrichment analysis by exposure category...")
  
  if (cor_df == "degs"){
    message("Enriching Differential Abundance Features Correlated with Exposure")
    correlation_results <- metadata(expOmicSet)$omics_exposure_deg_correlation
  } else if (cor_df == "factors"){
    message("Enriching Factor Features Correlated with Exposure")
    correlation_results <- metadata(expOmicSet)$omics_exposure_factor_correlation
  } else {
    stop("Invalid correlation dataframe specified. Use 'degs' or 'factors'.")
  }
  
  if (is.null(correlation_results) || nrow(correlation_results) == 0) {
    stop("No significant omics-exposure correlations found.")
  }
  
  unique_categories <- unique(correlation_results$category)
  enrichment_results <- list()
  
  for (cat in unique_categories) {
    message("Processing exposure category: ", cat)
    
    # **Extract assays linked to this exposure category**
    assay_names <- correlation_results |>
      filter(category == cat) |> 
      pull(assay_name) |> 
      unique()
    
    for (assay in assay_names) {
      message("Processing category: ", cat, " (Assay: ", assay, ")")
      
      # **Extract significant correlated features for this assay**
      sig_features <- correlation_results |> 
        filter(category == cat, 
               assay_name == assay,
               abs(correlation) > cor_threshold, 
               !!sym(pval_col) < pval_threshold) |> 
        pull(feature) |> 
        unique()
      
      if (length(sig_features) == 0) {
        message("No significant features for ", assay, " in category: ", cat, ". Skipping...")
        next
      }
      
      # **Dynamic Universe Setup (Skip for miRNAs)**
      tested_features <- if (!is.null(mirna_assays) && assay %in% mirna_assays) {
        NULL
      } else {
        rownames(experiments(expOmicSet)[[assay]])
      }
      
      # **Handle Protein Assays (Convert UniProt to Gene Symbols)**
      if (!is.null(uniprot_assays) && assay %in% uniprot_assays) {
        message("Converting UniProt IDs to gene symbols for: ", assay)
        sig_features <- as.character(sig_features)  # Ensure it's a character vector
        
        if (length(sig_features) > 0) {
          id_map <- .convert_uniprot_to_symbol(sig_features)
          sig_features <- id_map$hgnc_symbol[!is.na(id_map$hgnc_symbol)]
        }
        
        if (length(sig_features) == 0) {
          message("No valid gene symbols found after UniProt conversion for ", assay, ". Skipping...")
          next
        }
        
        # Convert the tested background features too
        if (!is.null(tested_features)) {
          tested_features <- .convert_uniprot_to_symbol(tested_features)$hgnc_symbol
        }
      }
      
      # **Handle miRNA Assays (Retrieve Target Genes)**
      if (!is.null(mirna_assays) && assay %in% mirna_assays) {
        message("Retrieving target genes for miRNAs in: ", assay)
        sig_features <- .get_mirna_targets(sig_features)
        
        if (length(sig_features) == 0) {
          message("No validated miRNA targets found for ", assay, ". Skipping...")
          next
        }
      }
      
      # **Ensure valid genes for enrichment**
      if (length(sig_features) == 0) {
        message("No valid genes for enrichment in assay: ", assay, ". Skipping...")
        next
      }
      
      # **Run GO Enrichment**
      message("Running GO enrichment for category: ", cat, " (Assay: ", assay, ", Ontology: ", ontology, ")")
      enrich_res <- tryCatch(
        enrichGO(
          gene = sig_features, 
          OrgDb = orgdb, 
          keyType = keytype,
          ont = ontology, 
          universe = if (is.null(tested_features) || !universe_background) NULL else tested_features,
          pAdjustMethod = p_adjust_method, 
          pvalueCutoff = pvalue_cutoff, 
          qvalueCutoff = qvalue_cutoff
        ),
        error = function(e) {
          message("GO enrichment failed for ", assay, " in category ", cat, ": ", e$message)
          return(NULL)
        }
      )
      
      # **Store results if enrichment was successful**
      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        enrich_df <- enrich_res@result |> mutate(category = cat, assay_name = assay)
        enrichment_results[[paste(cat, assay, sep = "_")]] <- enrich_df
      } else {
        message("No significant enrichment results for assay: ", assay, " in category: ", cat, ".")
      }
    }
  }
  
  # **Concatenate Results Across Categories & Assays**
  if (length(enrichment_results) > 0) {
    final_results <- bind_rows(enrichment_results)
    metadata(expOmicSet)$exposure_category_enrichment <- final_results
    metadata(expOmicSet)$exposure_category_enrichment_summary <- final_results |> 
      .summarize_exposure_enrichment()
  } else {
    message("No significant enrichment results found across all categories.")
    metadata(expOmicSet)$exposure_category_enrichment <- NULL
  }
  
  message("Enrichment analysis by exposure category completed.")
  return(expOmicSet)
}