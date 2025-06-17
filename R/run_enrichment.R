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
#' @return If `action="add"`, returns the updated `expomicset`.
#' If `action="get"`, returns a `data.frame` with enrichment results.
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
    action = "add"
) {
  # Run enrichment based on geneset type
  enrich_res <- switch(
    geneset,
    "deg" = .da_functional_enrichment(
      expomicset = expomicset,
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
    ),
    "deg_exp_cor" = .da_exposure_functional_enrichment(
      expomicset = expomicset,
      #feature_type = "degs",
      feature_col = feature_col,
      mirna_assays = mirna_assays,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      fun = fun,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont
    ),
    "factor_exp_cor" = .correlation_functional_enrichment(
      expomicset = expomicset,
      feature_type = "factors",
      feature_col = feature_col,
      mirna_assays = mirna_assays,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      fun = fun,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont
    ),
    "factor" = .factor_functional_enrichment(
      expomicset = expomicset,
      mirna_assays = mirna_assays,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont
    ),
    stop("Invalid geneset. Choose from 'deg', 'deg_exp_cor', 'factor_exp_cor', or 'factor'.")
  )

  if (is.null(enrich_res) || nrow(enrich_res) == 0) {
    warning("No enrichment results found for geneset '", geneset, "'. Skipping GO term clustering.")
    if (action == "add") {
      MultiAssayExperiment::metadata(expomicset)$enrichment[[geneset]] <- enrich_res
      return(expomicset)
    } else {
      return(enrich_res)
    }
  }

  message("Determining Number of GO Term Clusters...")

  # Deduplicate GO terms by Description
  go_term_list <- enrich_res |>
    dplyr::mutate(gene_list = stringr::str_split(geneID, "/")) |>
    tidyr::unnest(gene_list) |>
    dplyr::mutate(gene_list = stringr::str_trim(gene_list)) |>
    dplyr::distinct(Description, gene_list) |>
    dplyr::group_by(Description) |>
    dplyr::summarise(genes = list(unique(gene_list)), .groups = "drop")

  # Compute Jaccard matrix and cluster
  go_groups <- go_term_list$genes |>
    rlang::set_names(go_term_list$Description) |>
    .get_pairwise_overlaps() |>
    dplyr::select(source, target, jaccard) |>
    tidyr::pivot_wider(names_from = "target", values_from = "jaccard") |>
    tibble::column_to_rownames("source") |>
    as.matrix() |>
    .cluster_mat(clustering_approach = clustering_approach) |>
    (\(x) data.frame(Description = names(x), go_group = paste0("Group_", as.numeric(x))))()

  if (is.null(go_groups) || nrow(go_groups) == 0) {
    warning("GO clustering failed — no clusters found.")
    enrich_res$go_group <- NA_character_
  } else {
    enrich_res <- dplyr::inner_join(enrich_res, go_groups, by = "Description")
  }

  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$enrichment[[geneset]] <- enrich_res

    step_record <- list(
      run_enrichment = list(
        timestamp = Sys.time(),
        params = list(
          geneset = geneset,
          enrichment_fun = fun,
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff,
          clustering_approach = clustering_approach
        ),
        notes = paste0("Performed enrichment on ", geneset, " features.")
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)
  } else if (action == "get") {
    return(enrich_res)
  } else {
    stop("Invalid action. Choose from 'add' or 'get'.")
  }
}
