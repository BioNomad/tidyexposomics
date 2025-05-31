#' Perform Functional Enrichment Analysis
#'
#' Runs functional enrichment analysis on differentially abundant features or correlated features in a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param geneset A character string specifying the gene set to analyze. Options: `"deg"`, `"deg_exp_cor"`, `"factor_exp_cor"`, or `"factor"`.
#' @param feature_col A character string specifying the column for feature names. Default is `"feature"`.
#' @param mirna_assays A character vector specifying miRNA assays. Default is `NULL`.
#' @param pval_col A character string specifying the column for p-values. Default is `"adj.P.Val"`.
#' @param pval_threshold A numeric threshold for adjusted p-values. Default is `0.05`.
#' @param logfc_col A character string specifying the column for log-fold changes. Default is `"logFC"`.
#' @param logfc_threshold A numeric threshold for log-fold changes. Default is `log2(1)`.
#' @param pAdjustMethod A character string specifying the p-value adjustment method. Default is `"fdr"`.
#' @param pvalueCutoff A numeric threshold for significance in enrichment analysis. Default is `0.05`.
#' @param qvalueCutoff A numeric threshold for false discovery rate control. Default is `0.1`.
#' @param fun A character string specifying the enrichment function (`"enrichGO"`, `"enrichKEGG"`, etc.). Default is `"enrichGO"`.
#' @param OrgDb A character string specifying the organism database. Default is `"org.Hs.eg.db"`.
#' @param keyType A character string specifying the gene identifier type. Default is `"SYMBOL"`.
#' @param ont A character string specifying the ontology for Gene Ontology analysis. Default is `"BP"` (Biological Process).
#' @param clustering_approach A character string specifying the clustering approach for GO term grouping. Default is `"diana"`.
#' @param action A character string specifying `"add"` (store results in metadata) or `"get"` (return results). Default is `"add"`.
#'
#' @details
#' This function:
#' - Performs functional enrichment analysis on `geneset`-specific features.
#' - Supports enrichment for differentially abundant (`"deg"`) or correlated (`"deg_exp_cor"`, `"factor_exp_cor"`, `"factor"`) features.
#' - Uses `.get_pairwise_overlaps()` to cluster GO terms based on shared genes.
#' - Uses `.cluster_mat()` to determine GO term clusters using the selected `clustering_approach`.
#' - Stores results in `metadata(expomicset)$functional_enrichment[[geneset]]` when `action="add"`.
#'
#' @return If `action="add"`, returns the updated `expomicset`.
#' If `action="get"`, returns a `data.frame` with enrichment results, including:
#' \item{Description}{GO term description.}
#' \item{geneID}{Gene set associated with the GO term.}
#' \item{go_group}{Assigned GO term cluster.}
#'
#' @examples
#' \dontrun{
#' expom <- run_enrichment(
#'   expomicset = expom,
#'   geneset = "deg",
#'   pvalueCutoff = 0.01,
#'   clustering_approach = "gap"
#' )
#' }
#'
#' @export
run_enrichment <- function(
    expomicset,
    geneset,
    feature_col = "feature",
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
    clustering_approach = "diana",
    action="add"
){

  if(geneset == "deg"){
    enrich_res <- expomicset |>
      .da_functional_enrichment(
        mirna_assays = mirna_assays,
        pval_col = pval_col,
        pval_threshold = pval_threshold,
        logfc_col = logfc_col,
        logfc_threshold = logfc_threshold,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else if (geneset == "deg_exp_cor"){
    enrich_res <- expomicset |>
      .da_exposure_functional_enrichment(
        feature_col = feature_col,
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
    enrich_res <- expomicset |>
      .factor_exposure_functional_enrichment(
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
    enrich_res <- expomicset |>
      .factor_functional_enrichment(
        mirna_assays = mirna_assays,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont
      )
  } else {
    stop("Invalid geneset. Choose from 'deg', 'deg_cor', 'factor_exp_cor', or 'factor'.")
  }

  # Exit early if no enrichment was found
  if (is.null(enrich_res) || nrow(enrich_res) == 0) {
    warning("No enrichment results found for geneset '", geneset, "'. Skipping GO term clustering.")
    if (action == "add") {
      MultiAssayExperiment::metadata(expomicset)$functional_enrichment[[geneset]] <- enrich_res
      return(expomicset)
    } else {
      return(enrich_res)
    }
  }

  message("Determining Number of GO Term Clusters...")

  # determine go groups
  go_groups <- enrich_res |>
    (\(df) split(df,df$Description) )() |>
    purrr::map(~.x |>
          dplyr::pull(geneID) |>
          paste(collapse ="/") |>
          stringr::str_split("/") |>
          unlist() |>
          stringr::str_trim() |>
          unique()) |>
    .get_pairwise_overlaps() |>
    dplyr::select(source,target,jaccard) |>
    tidyr::pivot_wider(names_from = "target",
                values_from = "jaccard") |>
    tibble::column_to_rownames("source") |>
    as.matrix() |>
    .cluster_mat(clustering_approach = clustering_approach) |>
    (\(x) df=data.frame(
      Description=names(x),
      go_group=as.numeric(x)))()

  if (is.null(go_groups) || nrow(go_groups) == 0 || !"Description" %in% colnames(go_groups)) {
    warning("GO clustering failed — no Description terms or clusters found.")
    enrich_res$go_group <- NA_character_
    if (action == "add") {
      MultiAssayExperiment::metadata(expomicset)$functional_enrichment[[geneset]] <- enrich_res

      # Add analysis steps taken to metadata
      MultiAssayExperiment::metadata(expomicset)$steps <- c(
        MultiAssayExperiment::metadata(expomicset)$steps,
        "run_enrichment"
      )
      return(expomicset)
    } else {
      return(enrich_res)
    }
  }


  enrich_res <- enrich_res |>
    dplyr::inner_join(
      go_groups,
      by = "Description") |>
    dplyr::mutate(go_group=paste("Group", go_group, sep="_"))

  if(action=="add"){
    MultiAssayExperiment::metadata(expomicset)$functional_enrichment[[geneset]] <- enrich_res

    return(expomicset)
  }else if (action=="get"){
    return(enrich_res)
  }else{
    stop("Invalid action. Choose from 'add' or 'get'.")
  }
}
