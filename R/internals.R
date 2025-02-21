# --- Update Assay ColData ------
.update_assay_colData <- function(expOmicSet, exp_name) {
  require(MultiAssayExperiment)
  require(tidyverse)
  
  # Retrieve the assay
  assay <- experiments(expOmicSet)[[exp_name]]
  
  # Extract colData for the assay's samples
  assay_samples <- colnames(assay)
  global_coldata <- as.data.frame(colData(expOmicSet))
  coldata <- global_coldata[rownames(global_coldata) %in% assay_samples, , drop = FALSE]
  
  # Ensure the sample order matches
  coldata <- coldata[match(assay_samples, rownames(coldata)), , drop = FALSE]
  
  # Add a check to ensure the order is correct
  if (!identical(rownames(coldata), assay_samples)) {
    stop("Sample order mismatch detected in assay: ", exp_name, 
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
          feature = rownames(contrast_results),
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
        feature = rownames(contrast_results),
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
    group_by(feature, exp_name) |>
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
  require(biomaRt)
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
    dplyr::select(source_genesymbol,target_genesymbol)
  return(mirna_targets)
}
# --- Differential Abundance Exposure Functional Enrichment ------------
.da_exposure_functional_enrichment <- function(
    expOmicSet,
    proteomics_assays = NULL, 
    mirna_assays = NULL,
    pval_col = "adj.P.Val",
    pval_threshold = 0.05,
    logfc_col = "logFC",
    logfc_threshold = log2(1),
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {
  require(clusterProfiler)
  require(biomaRt)
  require(tidyverse)
  require(org.Hs.eg.db)
  
  
  if (!"differential_abundance" %in% names(expOmicSet@metadata)) {
    stop("Please run `run_differential_abundance() first.`")
  }
  
  if (!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)) {
    stop("Please run `correlate_exposures_with_degs() first.`")
  }
  
  da_res <- expOmicSet@metadata$differential_abundance
  da_cor_res <- expOmicSet@metadata$omics_exposure_deg_correlation
  
  da_res_cor_merged <- da_res |>
    filter(!!sym(pval_col) < pval_threshold) |>
    filter(abs(!!sym(logfc_col)) > logfc_threshold) |>
    mutate(direction=ifelse(logFC>0,"up","down")) |>
    dplyr::select(feature,
                  exp_name,
                  direction) |>
    inner_join(da_cor_res,
               by=c("feature",
                    "exp_name"))
  
  
  mirna_df <- da_res_cor_merged |>
    filter(exp_name %in% mirna_assays) |>
    (\(df) split(df,df$exp_name))() |>
    map(~ .get_mirna_targets(.x |> 
                               pull(feature)) |>
          mutate(direction="down")) |>
    bind_rows(.id="exp_name") |>
    inner_join(da_res_cor_merged |>
                 dplyr::select(-direction),
               by=c("source_genesymbol"="feature",
                    "exp_name")) |>
    dplyr::select(-source_genesymbol) |>
    dplyr::rename(feature=target_genesymbol)
  
  
  
  prot_df <- .convert_uniprot_to_symbol(da_res_cor_merged |>
                                          filter(exp_name %in% proteomics_assays) |>
                                          pull(feature)) |>
    inner_join(da_res_cor_merged,
               by=c("uniprot_gn_id"="feature")) |>
    dplyr::select(-uniprot_gn_id) |>
    dplyr::rename(feature=hgnc_symbol)
  
  
  da_res_cor_merged_mapped <- da_res_cor_merged |>
    filter(!exp_name %in% c(mirna_assays,proteomics_assays)) |>
    bind_rows(mirna_df,prot_df)
  
  
  universe_per_assay <- as.list(unique(da_res_cor_merged$exp_name)) |>
    map( ~
           {df <- data.frame(all_features=
                               experiments(expOmicSet)[[.x]] |>
                               rownames(),
                             exp_name=.x)
           
           if (.x %in% mirna_assays) {
             all_targets <- .get_mirna_targets(df$all_features) |>
               pull(target_genesymbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           } else if (.x %in% proteomics_assays){
             all_targets <- .convert_uniprot_to_symbol(df$all_features) |>
               pull(hgnc_symbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           }else{
             df <- df
           }
           return(df)}
    ) |>
    bind_rows()
  
  
  enrich_res <- da_res_cor_merged_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- compareCluster(
        feature~direction+category,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          filter(exp_name==unique(.x$exp_name)) |>
          pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })
  
  
  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")
  
  return(enrich_df)
}

# --- Differential Abundance Functional Enrichment ------------
.da_functional_enrichment <- function(
    expOmicSet,
    proteomics_assays = NULL, 
    mirna_assays = NULL,
    pval_col = "adj.P.Val",
    pval_threshold = 0.05,
    logfc_col = "logFC",
    logfc_threshold = log2(1),
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {
  require(clusterProfiler)
  require(biomaRt)
  require(tidyverse)
  require(org.Hs.eg.db)
  
  
  if (!"differential_abundance" %in% names(expOmicSet@metadata)) {
    stop("Please run `run_differential_abundance() first.`")
  }
  
  da_res <- expOmicSet@metadata$differential_abundance |>
    filter(!!sym(pval_col) < pval_threshold) |>
    filter(abs(!!sym(logfc_col)) > logfc_threshold) |>
    mutate(direction=ifelse(logFC>0,"up","down")) |>
    dplyr::select(feature,
                  exp_name,
                  direction) 
  
  
  mirna_df <- da_res |>
    filter(exp_name %in% mirna_assays) |>
    (\(df) split(df,df$exp_name))() |>
    map(~ .get_mirna_targets(.x |> 
                               pull(feature)) |>
          mutate(direction="down")) |>
    bind_rows(.id="exp_name") |>
    inner_join(da_res |>
                 dplyr::select(-direction),
               by=c("source_genesymbol"="feature",
                    "exp_name")) |>
    dplyr::select(-source_genesymbol) |>
    dplyr::rename(feature=target_genesymbol)
  
  
  
  prot_df <- .convert_uniprot_to_symbol(
    da_res |>
      filter(exp_name %in% proteomics_assays) |>
      pull(feature)) |>
    inner_join(da_res,
               by=c("uniprot_gn_id"="feature")) |>
    dplyr::select(-uniprot_gn_id) |>
    dplyr::rename(feature=hgnc_symbol)
  
  
  da_res_mapped <- da_res |>
    filter(!exp_name %in% c(mirna_assays,proteomics_assays)) |>
    bind_rows(mirna_df,prot_df)
  
  
  universe_per_assay <- as.list(unique(da_res$exp_name)) |>
    map( ~
           {df <- data.frame(all_features=
                               experiments(expOmicSet)[[.x]] |>
                               rownames(),
                             exp_name=.x)
           
           if (.x %in% mirna_assays) {
             all_targets <- .get_mirna_targets(df$all_features) |>
               pull(target_genesymbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           } else if (.x %in% proteomics_assays){
             all_targets <- .convert_uniprot_to_symbol(df$all_features) |>
               pull(hgnc_symbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           }else{
             df <- df
           }
           return(df)}
    ) |>
    bind_rows()
  
  
  enrich_res <- da_res_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- compareCluster(
        feature~direction,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          filter(exp_name==unique(.x$exp_name)) |>
          pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })
  
  
  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")
  
  return(enrich_df)
}

# --- Factor Feature Exposure Functional Enrichment ----------
.factor_exposure_functional_enrichment <- function(
    expOmicSet,
    proteomics_assays = NULL, 
    mirna_assays = NULL,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {
  require(clusterProfiler)
  require(biomaRt)
  require(tidyverse)
  require(org.Hs.eg.db)
  
  
  if (!"top_factor_features" %in% names(expOmicSet@metadata)) {
    stop("Please run `extract_top_factor_features() first.`")
  }
  
  if (!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)) {
    stop("Please run `correlate_exposures_with_factors() first.`")
  }
  
  factor_res <- expOmicSet@metadata$top_factor_features
  factor_cor_res <- expOmicSet@metadata$omics_exposure_factor_correlation
  
  factor_res_cor_merged <- factor_res |>
    dplyr::select(feature,
                  exp_name) |>
    inner_join(factor_cor_res,
               by=c("feature",
                    "exp_name"))
  
  
  mirna_df <- factor_res_cor_merged |>
    filter(exp_name %in% mirna_assays) |>
    (\(df) split(df,df$exp_name))() |>
    map(~ .get_mirna_targets(.x |> 
                               pull(feature))) |>
    bind_rows(.id="exp_name") |>
    inner_join(factor_res_cor_merged,
               by=c("source_genesymbol"="feature",
                    "exp_name")) |>
    dplyr::select(-source_genesymbol) |>
    dplyr::rename(feature=target_genesymbol)
  
  
  
  prot_df <- .convert_uniprot_to_symbol(factor_res_cor_merged |>
                                          filter(exp_name %in% proteomics_assays) |>
                                          pull(feature)) |>
    inner_join(factor_res_cor_merged,
               by=c("uniprot_gn_id"="feature")) |>
    dplyr::select(-uniprot_gn_id) |>
    dplyr::rename(feature=hgnc_symbol)
  
  
  factor_res_cor_merged_mapped <- factor_res_cor_merged |>
    filter(!exp_name %in% c(mirna_assays,proteomics_assays)) |>
    bind_rows(mirna_df,prot_df)
  
  
  universe_per_assay <- as.list(unique(factor_res_cor_merged$exp_name)) |>
    map( ~
           {df <- data.frame(all_features=
                               experiments(expOmicSet)[[.x]] |>
                               rownames(),
                             exp_name=.x)
           
           if (.x %in% mirna_assays) {
             all_targets <- .get_mirna_targets(df$all_features) |>
               pull(target_genesymbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           } else if (.x %in% proteomics_assays){
             all_targets <- .convert_uniprot_to_symbol(df$all_features) |>
               pull(hgnc_symbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           }else{
             df <- df
           }
           return(df)}
    ) |>
    bind_rows()
  
  
  enrich_res <- factor_res_cor_merged_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- compareCluster(
        feature~category,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          filter(exp_name==unique(.x$exp_name)) |>
          pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })
  
  
  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")
  
  return(enrich_df)
}

# --- Factor Feature Functional Enrichment ------------
.factor_functional_enrichment <- function(
    expOmicSet,
    proteomics_assays = NULL, 
    mirna_assays = NULL,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {
  require(clusterProfiler)
  require(biomaRt)
  require(tidyverse)
  require(org.Hs.eg.db)
  
  
  if (!"top_factor_features" %in% names(expOmicSet@metadata)) {
    stop("Please run `extract_top_factor_features() first.`")
  }
  
  factor_res <- expOmicSet@metadata$top_factor_features |>
    dplyr::select(feature,
                  exp_name) 
  
  
  mirna_df <- factor_res |>
    filter(exp_name %in% mirna_assays) |>
    (\(df) split(df,df$exp_name))() |>
    map(~ .get_mirna_targets(.x |> 
                               pull(feature))) |>
    bind_rows(.id="exp_name") |>
    inner_join(factor_res,
               by=c("source_genesymbol"="feature",
                    "exp_name")) |>
    dplyr::select(-source_genesymbol) |>
    dplyr::rename(feature=target_genesymbol)
  
  
  
  prot_df <- .convert_uniprot_to_symbol(factor_res |>
                                          filter(exp_name %in% proteomics_assays) |>
                                          pull(feature)) |>
    inner_join(factor_res,
               by=c("uniprot_gn_id"="feature")) |>
    dplyr::select(-uniprot_gn_id) |>
    dplyr::rename(feature=hgnc_symbol)
  
  
  factor_res_mapped <- factor_res |>
    filter(!exp_name %in% c(mirna_assays,proteomics_assays)) |>
    bind_rows(mirna_df,prot_df)
  
  
  universe_per_assay <- as.list(unique(factor_res$exp_name)) |>
    map( ~
           {df <- data.frame(all_features=
                               experiments(expOmicSet)[[.x]] |>
                               rownames(),
                             exp_name=.x)
           
           if (.x %in% mirna_assays) {
             all_targets <- .get_mirna_targets(df$all_features) |>
               pull(target_genesymbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           } else if (.x %in% proteomics_assays){
             all_targets <- .convert_uniprot_to_symbol(df$all_features) |>
               pull(hgnc_symbol)
             
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           }else{
             df <- df
           }
           return(df)}
    ) |>
    bind_rows()
  
  
  enrich_res <- factor_res_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- enrichGO(
        gene = .x$feature,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          filter(exp_name==unique(.x$exp_name)) |>
          pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      
      )
    })
  
  
  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@result
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")
  
  return(enrich_df)
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
          mutate(exp_name = assay)
      } else {
        NULL
      }
    }) |> 
    bind_rows() |> 
    separate_rows(geneID, sep = "/") |> 
    group_by(ID, Description, exp_name) |> 
    reframe(
      genes = paste(unique(geneID), collapse = ", "), 
      num_genes = n_distinct(geneID),  # Count of unique genes per omic
      min_p_adj = min(p.adjust, na.rm = TRUE),
      max_p_adj = max(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    group_by(ID, Description) |> 
    mutate(n_omics = n_distinct(exp_name)) |>  # Calculate before using it
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
# --- Plot Pie Chart ----------------
.plot_circular_bar <- function(data, category_col, count_col) {
  data |> 
    mutate(
      fraction = !!sym(count_col) / sum(!!sym(count_col)),  # Compute percentages
      ymax = cumsum(fraction),  # Top of each rectangle
      ymin = c(0, head(ymax, n = -1)),  # Bottom of each rectangle
      labelPosition = (ymax+ymin)/2,  # Midpoint for labels
      label = paste0(!!sym(category_col),":", " ", !!sym(count_col))  # Label formatting
    ) |> 
    ggplot(aes(
      ymax = ymax,
      ymin = ymin,
      xmax = 4, 
      xmin = 3,
      fill = !!sym(category_col)
    )) +
    geom_rect(alpha=0.8) +
    geom_label_repel(
      x = 2, aes(y = labelPosition, 
                 label = label
                 #color = !!sym(category_col)
                 ), 
      size = 4,
      fill = "white",
      color="black") +
    scale_fill_cosmic()+
    coord_polar(theta = "y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none")
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
    dplyr::select(ID, Description, p.adjust, geneID, category, exp_name) |>
    separate_rows(geneID, sep = "/") |>  # Split gene lists into separate rows
    group_by(ID, Description, category, exp_name) |> 
    reframe(
      genes = paste(unique(geneID), collapse = ", "), 
      num_genes = n_distinct(geneID),  # Count of unique genes per omic-assay
      min_p_adj = min(p.adjust, na.rm = TRUE),
      max_p_adj = max(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    group_by(ID, Description, category) |> 
    mutate(n_omics = n_distinct(exp_name)) |>  # Count how many omics report the GO term
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
# --- Pairwise Overlaps -------------

.get_pairwise_overlaps <- function(sets) {
  # credit for most of the code:
  # https://blog.jdblischak.com/posts/pairwise-overlaps/
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x) length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }
  
  n_sets <- length(sets)
  set_names <- names(sets)
  n_overlaps <- choose(n = n_sets, k = 2)
  
  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_num_shared <- integer(length = n_overlaps)
  vec_overlap <- numeric(length = n_overlaps)
  vec_jaccard <- numeric(length = n_overlaps)
  vec_shared_terms <- character(length = n_overlaps)
  overlaps_index <- 1
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]
      
      shared_terms <- paste(Reduce(intersect,list(set1,set2)),collapse = ",")
      
      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                    fromLast = FALSE, nmax = NA))
      num_shared <- length(set_intersect)
      overlap <- num_shared / min(length(set1), length(set2))
      jaccard <- num_shared / length(set_union)
      
      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_num_shared[overlaps_index] <- num_shared
      vec_overlap[overlaps_index] <- overlap
      vec_jaccard[overlaps_index] <- jaccard
      vec_shared_terms[overlaps_index] <- shared_terms
      
      overlaps_index <- overlaps_index + 1
    }
  }
  
  result <- data.frame(source = vec_name1,
                       target = vec_name2,
                       num_shared = vec_num_shared,
                       overlap = vec_overlap,
                       jaccard = vec_jaccard,
                       shared_terms = vec_shared_terms,
                       stringsAsFactors = FALSE)
  return(result)
}
# --- Cluster Matrix ------------------
.cluster_mat <- function(data_matrix, dist_method = NULL, cluster_method = "ward.D", clustering_approach = "gap") {
  require(tidyverse)
  require(ComplexHeatmap)
  require(circlize)
  require(cluster)  # For silhouette scores
  require(vegan)  # For Gower's distance
  require(FactoMineR)  # For handling categorical data
  require(factoextra)
  require(dynamicTreeCut)
  require(densityClust)
  
  # Determine appropriate distance metric
  if (is.null(dist_method)) {
    if (all(sapply(data_matrix, is.numeric))) {
      dist_method <- "euclidean"  # Continuous data
    } else {
      dist_method <- "gower"  # Mixed data types
    }
  }
  
  # Compute distance matrix
  sample_dist <- if (dist_method == "gower") {
    daisy(data_matrix, metric = "gower")
  } else {
    dist(data_matrix, method = dist_method)
  }
  
  # Function to determine optimal k based on the selected clustering approach
  determine_k <- function(dist_matrix, cluster_method) {
    if (clustering_approach == "diana") {
      sample_cluster <- diana(as.dist(dist_matrix))
      height_diffs <- diff(sample_cluster$height)
      cutoff_index <- which.max(height_diffs)
      return(length(sample_cluster$height) - cutoff_index)
      
    } else if (clustering_approach == "gap") {
      gap_stat <- clusGap(as.matrix(dist_matrix), FUN = hcut, K.max = 20, B = 50)
      return(maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"]))
      
    } else if (clustering_approach == "elbow") {
      fviz_nbclust(as.matrix(dist_matrix), FUN = hcut, method = "wss")
      return(3)  # Adjust manually if needed
      
    } else if (clustering_approach == "dynamic") {
      sample_cluster <- hclust(as.dist(dist_matrix), method = cluster_method)
      cut_clusters <- cutreeDynamic(dendro = as.dendrogram(sample_cluster), distM = as.matrix(dist_matrix), deepSplit = 2)
      return(length(unique(cut_clusters)))
      
    } else if (clustering_approach == "density") {
      dclust <- densityClust(dist_matrix, gaussian = TRUE)
      dclust <- findClusters(dclust, rho = 0.3, delta = 0.5)
      return(max(dclust$clusters))
      
    } else {
      stop("Invalid clustering approach selected.")
    }
  }
  
  k_samples <- determine_k(sample_dist, cluster_method)
  
  # Perform hierarchical clustering only if needed
  sample_cluster <- hclust(as.dist(sample_dist), method = cluster_method)
  
  # Cut dendrograms using optimal k values
  sample_groups <- cutree(sample_cluster, k = k_samples)
  
  message("Optimal number of clusters for samples: ", k_samples)
  
  return(sample_groups)
}
# --- Test DA Functional Enrichment ------------

# --- Next Function --------