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
    ont = "BP"
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
  
  # create a list with the results
  enrich_res_lst <- list(enrich_res)
  names(enrich_res_lst) <- geneset
  
  expOmicSet@metadata$functional_enrichment <- enrich_res_lst
  
  return(expOmicSet)
}
