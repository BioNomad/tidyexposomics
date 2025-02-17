# --- Update Assay ColData ------
.update_assay_colData <- function(expOmicSet, assay_name) {
  require(MultiAssayExperiment)
  require(tidyverse)
  
  # Retrieve the assay
  assay <- experiments(expOmicSet)[[assay_name]]
  
  # Extract colData for the assay's samples
  assay_samples <- colnames(assay)
  global_coldata <- as.data.frame(colData(expOmicSet))
  coldata <- global_coldata[rownames(global_coldata) %in% assay_samples, , drop = FALSE]
  
  # Ensure the sample order matches
  coldata <- coldata[match(assay_samples, rownames(coldata)), , drop = FALSE]
  
  # Add a check to ensure the order is correct
  if (!identical(rownames(coldata), assay_samples)) {
    stop("Sample order mismatch detected in assay: ", assay_name, 
         "\nEnsure the samples in colData are aligned with the assay samples.")
  }
  
  # Update colData in the assay
  colData(assay) <- DataFrame(coldata)
  
  return(assay)
}


# --- Run Differential Abundance Analysis ------

.run_se_differential_abundance <- function(
    exp,  # SummarizedExperiment object
    formula,  # Model formula
    abundance_col = "counts",
    method = "limma_voom",
    scaling_method = "none",
    min_counts = 10,
    min_proportion = 0.3,
    contrasts = NULL
) {
  require(tidybulk)
  require(tidyverse)
  
  # --- Run differential expression analysis -----
  if (!is.null(contrasts)) {
    res_list <- list()
    for (contrast in contrasts) {
      contrast_results <- exp |>
        identify_abundant(minimum_counts = min_counts,
                          minimum_proportion = min_proportion) |>
        test_differential_abundance(
          formula,
          .abundance = !!sym(abundance_col),
          method = method,
          contrasts = contrast,
          scaling_method = scaling_method
        )
      
      # Extract results
      res <- as.data.frame(elementMetadata(contrast_results))
      colnames(res) <- gsub("__.*", "", colnames(res))
      
      # Add metadata
      res <- res |> 
        mutate(
          molecular_feature = rownames(contrast_results),
          contrast = contrast,
          method = method,
          scaling = scaling_method,
          min_counts = min_counts,
          min_proportion = min_proportion
        )
      
      res_list[[contrast]] <- res
    }
    return(bind_rows(res_list))
  } else {
    contrast_results <- exp |>
      identify_abundant(minimum_counts = min_counts,
                        minimum_proportion = min_proportion) |>
      test_differential_abundance(
        formula,
        .abundance = !!sym(abundance_col),
        method = method,
        scaling_method = scaling_method
      )
    
    # Extract results
    res <- as.data.frame(elementMetadata(contrast_results))
    colnames(res) <- gsub("__.*", "", colnames(res))
    
    # Add metadata
    res <- res |> 
      mutate(
        molecular_feature = rownames(contrast_results),
        contrast = all.vars(formula)[1],
        method = method,
        scaling = scaling_method,
        min_counts = min_counts,
        min_proportion = min_proportion
      )
    
    return(res)
  }
}

# --- Calculate Feature Stability -------
.calculate_feature_stability <- function(sensitivity_df, 
                                         pval_col = "adj.P.Val",
                                         logfc_col = "logFC",
                                         pval_threshold = 0.05, 
                                         logFC_threshold = log2(1.5)) {
  # Extract sensitivity analysis results
  sensitivity_df <- sensitivity_df
  
  message("Computing feature stability across sensitivity conditions...")
  
  # Identify significant features in each test setting
  feature_stability_df <- sensitivity_df |>
    filter(!!sym(pval_col) < pval_threshold,
           abs(!!sym(logfc_col)) > logFC_threshold) |>
    group_by(molecular_feature, exp_name) |>
    dplyr::summarize(stability_score = n(), .groups = "drop") |>
    arrange(desc(stability_score))
  
  # Store results in metadata
  # metadata(expOmicSet)$feature_stability_df <- feature_stability_df
  
  message("Feature stability analysis completed.")
  return(feature_stability_df)
}

# --- Correlate Se with colData --------

.correlate_se_with_coldata <- function(
    se,  # SummarizedExperiment object
    exposure_cols,  # Character vector of colData columns for correlation
    correlation_method = "spearman",
    correlation_cutoff = 0.0,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05
) {
  require(tidyverse)
  require(Hmisc)
  require(reshape2)
  
  message("Performing correlation analysis on summarized experiment...")
  
  # Ensure colData has the specified exposures
  exposure_data <- colData(se) |> as.data.frame()
  numeric_exposures <- intersect(colnames(exposure_data), exposure_cols)
  
  if (length(numeric_exposures) == 0) {
    stop("No valid exposure variables found in colData.")
  }
  
  # Extract assay data
  assay_data <- assays(se)[[1]] |> t() |> as.data.frame()
  
  # Ensure colData and assay samples match
  common_samples <- intersect(rownames(assay_data), rownames(exposure_data))
  
  if (length(common_samples) == 0) {
    stop("No common samples between assay and colData.")
  }
  
  # Subset both data frames to common samples
  assay_data <- assay_data[common_samples, , drop = FALSE]
  exposure_data <- exposure_data[common_samples, numeric_exposures, drop = FALSE]
  
  # Merge into a single dataframe for correlation
  merged_data <- assay_data |> 
    rownames_to_column("id") |> 
    inner_join(exposure_data |> rownames_to_column("id"), by = "id") |> 
    column_to_rownames("id")
  
  # Perform correlation analysis
  message("Running Spearman correlation analysis...")
  correlation_matrix <- merged_data |> as.matrix() |> Hmisc::rcorr(type = correlation_method)
  
  # Convert correlation and p-values to tidy format
  correlation_df <- correlation_matrix$r |> 
    as.data.frame() |> 
    rownames_to_column("feature") |> 
    melt(id.vars = "feature") |> 
    `colnames<-`(c("feature", "exposure", "correlation"))
  
  pvalue_df <- correlation_matrix$P |> 
    as.data.frame() |> 
    rownames_to_column("feature") |> 
    melt(id.vars = "feature") |> 
    `colnames<-`(c("feature", "exposure", "p.value"))
  
  # Merge correlation results
  correlation_results <- inner_join(correlation_df,
                                    pvalue_df, 
                                    by = c("feature", "exposure")) |> 
    filter(!(abs(correlation) > 1)) |>  # Remove perfect correlations
    filter(abs(correlation) > correlation_cutoff) |>  # Apply correlation threshold
    mutate(FDR = p.adjust(p.value, method = "fdr")) |>  # Adjust p-values
    filter(!!sym(cor_pval_column) < pval_cutoff) |>  # Apply p-value cutoff
    arrange(desc(abs(correlation))) |>  # Sort by strongest correlations 
    filter(feature %in% rownames(se)) |>  # Filter to features in SummarizedExperiment
    filter(exposure %in% numeric_exposures)  # Filter to valid exposures
  
  message("Correlation analysis completed.")
  
  return(correlation_results)
}
# --- Convert Uniprot to Symbol ------

.convert_uniprot_to_symbol <- function(uniprot_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_map <- getBM(attributes = c("uniprot_gn_id", "hgnc_symbol"),
                    filters = "uniprot_gn_id",
                    values = uniprot_ids,
                    mart = mart)
  return(gene_map)
}

# --- Get miRNA Targets ----------

.get_mirna_targets <- function(mirnas) {
  # approach using multimir
  # targets <- get_multimir(org = "hsa", mirna = mirnas, table = "validated")
  # return(targets@data$target_symbol)
  
  # pull omnipath database
  require(OmnipathR)
  mirna_db <- OmnipathR::import_mirnatarget_interactions()
  mirna_targets <- mirna_db |>
    filter(source_genesymbol %in% mirnas) |>
    pull(target_genesymbol) |> 
    unique()
}
# --- Differential Abundance Functional Enrichment ------------
.da_functional_enrichment <- function(
    expOmicSet, 
    pval_col = "adj.P.Val", 
    logfc_col = "logFC", 
    pval_threshold = 0.05, 
    logFC_threshold = log2(1.5), 
    ontology = "BP",  # BP (biological process), MF (molecular function), CC (cellular component)
    mirna_assays = NULL,  # User-specified miRNA assay names
    uniprot_assays = NULL,  # User-specified protein assay names with UniProt IDs
    universe_background = TRUE,  # Whether to set universe background (disabled for miRNA)
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
  
  message("Starting differential abundance functional enrichment analysis...")
  
  da_results <- metadata(expOmicSet)$differential_abundance
  
  if (is.null(da_results) || nrow(da_results) == 0) {
    stop("No differential abundance results found in metadata.")
  }
  
  unique_assays <- unique(da_results$assay_name)
  enrichment_results <- list()
  
  for (assay in unique_assays) {
    message("Processing assay: ", assay)
    
    # Extract all features tested for dynamic universe (skip miRNA)
    tested_features <- if (!is.null(mirna_assays) && assay %in% mirna_assays) NULL else rownames(experiments(expOmicSet)[[assay]])
    
    # Extract significant DE features
    sig_features <- da_results |> 
      filter(assay_name == assay, 
             !!sym(pval_col) < pval_threshold, 
             abs(!!sym(logfc_col)) > logFC_threshold) |> 
      pull(molecular_feature) |> unique()
    
    if (length(sig_features) == 0) {
      message("No significant features for ", assay, ". Skipping...")
      next
    }
    
    # Handle Protein Assays (Convert UniProt to Gene Symbols)
    if (!is.null(uniprot_assays) && assay %in% uniprot_assays) {
      message("Converting UniProt IDs to gene symbols for: ", assay)
      id_map <- .convert_uniprot_to_symbol(sig_features)
      sig_features <- id_map$hgnc_symbol[!is.na(id_map$hgnc_symbol)]
      
      if (length(sig_features) == 0) {
        message("No valid gene symbols found after UniProt conversion for ", assay, ". Skipping...")
        next
      }
      
      # Convert tested features if universe_background is TRUE
      if (universe_background && !is.null(tested_features)) {
        tested_features <- .convert_uniprot_to_symbol(tested_features)$hgnc_symbol
      }
    }
    
    # Handle miRNA Assays (Get Target Genes)
    if (!is.null(mirna_assays) && assay %in% mirna_assays) {
      message("Retrieving target genes for miRNAs in: ", assay)
      sig_features <- .get_mirna_targets(sig_features)
      
      if (length(sig_features) == 0) {
        message("No validated miRNA targets found for ", assay, ". Skipping...")
        next
      }
    }
    
    # Ensure valid genes for enrichment
    if (length(sig_features) == 0) {
      message("No valid gene symbols found for ", assay, ". Skipping enrichment...")
      next
    }
    
    # Run GO Enrichment
    message("Running GO enrichment for ", assay, " (", ontology, ")")
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
        message("GO enrichment failed for ", assay, ": ", e$message)
        return(NULL)
      }
    )
    
    # Store results if enrichment was successful
    if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
      enrichment_results[[assay]] <- enrich_res@result
    } else {
      message("No significant enrichment results for assay: ", assay)
    }
  }
  
  # Store results in metadata
  metadata(expOmicSet)$da_enrichment_results <- enrichment_results
  message("Differential abundance functional enrichment analysis completed.")
  return(expOmicSet)
}


# --- Factor Feature Functional Enrichment ------------
.factor_functional_enrichment <- function(
    expOmicSet, 
    ontology = "BP",  # BP (biological process), MF (molecular function), CC (cellular component)
    min_loading = 0.0,  # Minimum absolute loading threshold for feature inclusion
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
  
  message("Starting functional enrichment analysis for factor features...")
  
  factor_features <- metadata(expOmicSet)$top_factor_features
  
  if (is.null(factor_features) || nrow(factor_features) == 0) {
    stop("No factor features found in metadata.")
  }
  
  unique_assays <- unique(factor_features$name)
  enrichment_results <- list()
  
  for (assay in unique_assays) {
    message("Processing assay: ", assay)
    
    # **Extract significant factor-contributing features**
    sig_features <- factor_features |> 
      filter(name == assay, abs(loading) > min_loading) |> 
      pull(features) |> unique()
    
    if (length(sig_features) == 0) {
      message("No significant factor-contributing features for ", assay, ". Skipping...")
      next
    }
    
    # **Extract all tested features for dynamic universe (skip miRNA)**
    tested_features <- if (!is.null(mirna_assays) && assay %in% mirna_assays) NULL else rownames(experiments(expOmicSet)[[assay]])
    
    # **Handle Protein Assays (Convert UniProt to Gene Symbols)**
    if (!is.null(uniprot_assays) && assay %in% uniprot_assays) {
      message("Converting UniProt IDs to gene symbols for: ", assay)
      id_map <- .convert_uniprot_to_symbol(sig_features)
      sig_features <- id_map$hgnc_symbol[!is.na(id_map$hgnc_symbol)]
      
      if (length(sig_features) == 0) {
        message("No valid gene symbols found after UniProt conversion for ", assay, ". Skipping...")
        next
      }
      
      # Convert full tested feature list to gene symbols if background is enabled
      if (universe_background && !is.null(tested_features)) {
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
    
    # **Ensure valid genes for enrichment**
    if (length(sig_features) == 0) {
      message("No valid gene symbols found for ", assay, ". Skipping enrichment...")
      next
    }
    
    # **Run GO Enrichment**
    message("Running GO enrichment for ", assay, " (", ontology, ")")
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
        message("GO enrichment failed for ", assay, ": ", e$message)
        return(NULL)
      }
    )
    
    # **Store results if enrichment was successful**
    if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
      enrichment_results[[assay]] <- enrich_res@result
    } else {
      message("No significant enrichment results for ", assay, ".")
    }
  }
  
  # **Store results in metadata**
  metadata(expOmicSet)$factor_feature_enrichment <- enrichment_results
  message("Factor-based functional enrichment analysis completed.")
  return(expOmicSet)
}

# --- Summarize Go Enrichment ---------------------
.summarize_go_enrichment <- function(enrichment_results, p_adjust_threshold = 0.05) {
  require(tidyverse)
  
  if (length(enrichment_results) == 0) {
    stop("No enrichment results available for summarization.")
  }
  
  message("Summarizing GO enrichment results...")
  
  # Convert enrichment list to a tidy dataframe
  enriched_df <- names(enrichment_results) |> 
    map(~ {
      assay <- .x
      if (!is.null(enrichment_results[[assay]]) && nrow(enrichment_results[[assay]]) > 0) {
        enrichment_results[[assay]] |>
          filter(p.adjust < p_adjust_threshold) |>  # Apply p.adjust threshold
          dplyr::select(ID, Description, p.adjust, geneID) |>
          mutate(assay_name = assay)
      } else {
        NULL
      }
    }) |> 
    bind_rows() |> 
    separate_rows(geneID, sep = "/") |> 
    group_by(ID, Description, assay_name) |> 
    reframe(
      genes = paste(unique(geneID), collapse = ", "), 
      num_genes = n_distinct(geneID),  # Count of unique genes per omic
      min_p_adj = min(p.adjust, na.rm = TRUE),
      max_p_adj = max(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    group_by(ID, Description) |> 
    mutate(n_omics = n_distinct(assay_name)) |>  # Calculate before using it
    mutate(
      all_genes = paste(unique(genes), collapse = ", "), 
      shared_genes = if (dplyr::first(n_omics) > 1) {
        shared_list <- Reduce(intersect, strsplit(unique(genes), ", "))
        if (length(shared_list) == 0) NA_character_ else paste(shared_list, collapse = ", ")
      } else {
        NA_character_
      },
      num_shared_genes = ifelse(is.na(shared_genes), 0, length(unique(unlist(strsplit(shared_genes, ", "))))),
      num_all_genes = n_distinct(unlist(strsplit(all_genes, ", "))),
      category = case_when(
        n_omics == 1 ~ "Unique",
        n_omics > 1 ~ "Shared"
      )
    ) |> ungroup()
  
  message("GO enrichment summary completed.")
  
  return(enriched_df)
}

# --- Summarize Exposure Enrichment Results -----------------------
.summarize_exposure_enrichment <- function(enrichment_df, p_adjust_threshold = 0.05) {
  require(tidyverse)
  
  if (nrow(enrichment_df) == 0) {
    stop("No enrichment results available for summarization.")
  }
  
  message("Summarizing GO enrichment results...")
  
  # Filter by p.adjust threshold
  enriched_df <- enrichment_df |>
    filter(p.adjust < p_adjust_threshold) |>
    dplyr::select(ID, Description, p.adjust, geneID, category, assay_name) |>
    separate_rows(geneID, sep = "/") |>  # Split gene lists into separate rows
    group_by(ID, Description, category, assay_name) |> 
    reframe(
      genes = paste(unique(geneID), collapse = ", "), 
      num_genes = n_distinct(geneID),  # Count of unique genes per omic-assay
      min_p_adj = min(p.adjust, na.rm = TRUE),
      max_p_adj = max(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    group_by(ID, Description, category) |> 
    mutate(n_omics = n_distinct(assay_name)) |>  # Count how many omics report the GO term
    mutate(
      all_genes = paste(unique(genes), collapse = ", "), 
      shared_genes = if (dplyr::first(n_omics) > 1) {
        shared_list <- Reduce(intersect, strsplit(unique(genes), ", "))
        if (length(shared_list) == 0) NA_character_ else paste(shared_list, collapse = ", ")
      } else {
        NA_character_
      },
      num_shared_genes = ifelse(is.na(shared_genes), 0, length(unique(unlist(strsplit(shared_genes, ", "))))),
      num_all_genes = n_distinct(unlist(strsplit(all_genes, ", "))),
      category_status = case_when(
        n_omics == 1 ~ "Unique",
        n_omics > 1 ~ "Shared"
      )
    ) |> ungroup()
  
  message("GO enrichment summary completed.")
  
  return(enriched_df)
}

# --- Next Function --------