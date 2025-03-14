#' Plot Exposures Associated with Features in GO Groups
#'
#' Visualizes exposures that are significantly associated with features in specified Gene Ontology (GO) groups.
#'
#' @param expomicset A `MultiAssayExperiment` object containing functional enrichment and exposure-feature correlation results.
#' @param go_groups A character vector specifying GO groups of interest. Use `"all"` to include all GO groups.
#' @param geneset A character string indicating the gene set to use. Options are `"deg_exp_cor"` for DEGs or `"factor_exp_cor"` for factors. Default is `"deg_exp_cor"`.
#' @param feature_col A character string specifying the feature column to match with correlated features. Default is `"feature"`.
#' @param top_n An integer specifying the number of top exposures to display per experimental assay. Default is `15`.
#'
#' @details
#' This function identifies exposures significantly correlated with features within specified GO groups. 
#' It extracts functional enrichment results from `metadata(expomicset)$functional_enrichment` and 
#' cross-references them with exposure-feature correlation results:
#'
#' - `"deg_exp_cor"`: Uses `metadata(expomicset)$omics_exposure_deg_correlation`
#' - `"factor_exp_cor"`: Uses `metadata(expomicset)$omics_exposure_factor_correlation`
#'
#' The resulting plot:
#' - Displays exposures grouped by experimental assay (`exp_name`).
#' - Uses bar heights to indicate the number of associated features.
#' - Colors bars based on exposure categories.
#'
#' @return A `ggplot` object displaying the top exposures associated with features in GO groups.
#'
#' @examples
#' \dontrun{
#' plot_go_group_exposures(expom, go_groups = c("GO:0006355", "GO:0045893"))
#' }
#'
#' @export
plot_go_group_exposures <- function(
    expomicset,
    go_groups,
    geneset = "deg_exp_cor",
    feature_col = "feature",
    top_n = 15
){
  require(ggplot2)
  
  # Check if the required metadata is present
  if(!"functional_enrichment" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_enrichment()` first.")
  }
  
  # Get the functional enrichment results
  enrich_res <- MultiAssayExperiment::metadata(expomicset)$functional_enrichment[[geneset]]
  
  # Filter the results to the GO groups of interest
  if(identical(go_groups,"all")){
    go_group_res <- enrich_res 
  }else{
    go_group_res <- enrich_res |> 
      dplyr::filter(go_group %in% go_groups)
  }
  
  # Get the genes in the GO groups
  go_group_genes <- go_group_res |> 
    dplyr::pull(geneID) |> 
    stringr::str_split("/") |> 
    unlist() |> 
    unique()
  
  # Get the correlation results
  if(geneset=="deg_exp_cor"){
    cor_res <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation |> 
      dplyr::inner_join(pivot_feature(expomicset),
                 by=c("feature"=".feature",
                      "exp_name"=".exp_name")) 
  }else if(geneset=="factor_exp_cor"){
    cor_res <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation |> 
      dplyr::inner_join(pivot_feature(expomicset),
                 by=c("feature"=".feature",
                      "exp_name"=".exp_name"))
  }
  
  # Filter the correlation results to the GO groups
  cor_res_go_group <- cor_res |> 
    dplyr::filter(!!sym(feature_col) %in% go_group_genes)
  
  # Plot the results
  cor_res_go_group |>
    dplyr::group_by(exp_name,category,exposure) |>
    dplyr::summarise(n=dplyr::n()) |>
    dplyr::ungroup() |> 
    dplyr::group_by(exp_name) |> 
    dplyr::arrange(desc(n)) |>
    dplyr::slice_head(n=top_n) |>
    dplyr::ungroup() |>
    ggplot(aes(
      x=n,
      y=tidytext::reorder_within(exposure,n,exp_name),
      fill=category,
      group=exp_name
    ))+
    geom_bar(stat="identity",alpha=0.7)+
    facet_wrap(exp_name~.,scales="free_y",ncol=1)+
    theme_bw()+
    theme(strip.text.x = element_text(angle=0,face = "bold.italic"),
          plot.subtitle = element_text(face="italic"),
          plot.title = element_text(face="bold.italic"))+
    ggsci::scale_fill_cosmic()+
    ggsci::scale_color_cosmic()+
    tidytext::scale_y_reordered()+
    labs(
      x="Number of Feature Associations",
      y="",
      fill="Category",
      title="Exposures Associated with Features in GO Groups",
      subtitle = paste("GO Groups:",paste(go_groups,collapse = ","))
    )

}

