run_functional_enrichment <- function(
    expOmicSet,
    geneset, # deg, factor, deg_exp_cor,factor_exp_cor
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
    ont = "BP",
    clustering_approach = "diana"
){
  
  if(geneset == "deg"){
    enrich_res <- expOmicSet |> 
      .da_functional_enrichment(
        proteomics_assays = proteomics_assays, 
        mirna_assays = mirna_assays,
        pval_col = pval_col,
        pval_threshold = pval_threshold,
        logfc_col = logfc_col,
        logfc_threshold = logfc_threshold,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else if (geneset == "deg_exp_cor"){
    enrich_res <- expOmicSet |> 
      .da_exposure_functional_enrichment(
        proteomics_assays = proteomics_assays, 
        mirna_assays = mirna_assays,
        pval_col = pval_col,
        pval_threshold = pval_threshold,
        logfc_col = logfc_col,
        logfc_threshold = logfc_threshold,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else if (geneset == "factor_exp_cor"){
    enrich_res <- expOmicSet |> 
      .factor_exposure_functional_enrichment(
        proteomics_assays = proteomics_assays, 
        mirna_assays = mirna_assays,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else if (geneset == "factor"){
    enrich_res <- expOmicSet |> 
      .factor_functional_enrichment(
        proteomics_assays = proteomics_assays, 
        mirna_assays = mirna_assays,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else {
    stop("Invalid geneset. Choose from 'deg', 'deg_cor', 'factor_exp_cor', or 'factor'.")
  }
  
  message("Determining Number of GO Term Clusters...")
  
  # determine go groups
  go_groups <- enrich_res |>
    (\(df) split(df,df$Description) )() |>
    map(~.x |>
          pull(geneID) |>
          paste(collapse ="/") |>
          str_split("/") |>
          unlist() |>
          str_trim() |>
          unique()) |>
    .get_pairwise_overlaps() |>
    dplyr::select(source,target,jaccard) |>
    pivot_wider(names_from = "target",
                values_from = "jaccard") |>
    column_to_rownames("source") |>
    as.matrix() |>
    .cluster_mat(clustering_approach = clustering_approach) |>
    (\(x) df=data.frame(
      Description=names(x),
      go_group=as.numeric(x)))()
  
  # create a list with the results
  enrich_res_lst <- list(
    enrich_res=enrich_res,
    go_groups=go_groups)
  
  expOmicSet@metadata$functional_enrichment[[geneset]] <- enrich_res_lst
  
  return(expOmicSet)
}
