#' Plot Dotplot of Functional Enrichment Results
#'
#' Generates a dotplot visualization of enriched gene ontology (GO) terms
#' across different experimental categories.
#'
#' @param expomicset A `MultiAssayExperiment` object containing functional enrichment results.
#' @param geneset A character string specifying the gene set name to extract enrichment results from.
#' @param top_n An integer specifying the number of top GO groups to display. Default is `5`.
#' @param n_per_group An integer specifying the number of top enriched terms to display within each GO group. Default is `5`.
#' @param go_groups A character vector specifying specific GO groups to include. If `NULL`, selects the top `top_n` groups based on enrichment score.
#'
#' @details
#' This function extracts functional enrichment results from `metadata(expomicset)$functional_enrichment`,
#' selects the most significant GO terms based on `-log10(p.adjust) * Count`,
#' and visualizes them in a dotplot.
#'
#' - The x-axis represents experimental categories.
#' - The y-axis represents enriched GO terms.
#' - Dot size indicates the number of genes in the GO term.
#' - Dot color represents statistical significance (`-log10(p.adjust)`).
#'
#' The function allows filtering by predefined `go_groups` or selecting the top `top_n` groups automatically.
#'
#' @return A `ggplot` object displaying a dotplot of enriched GO terms.
#'
#' @examples
#' \dontrun{
#' plot_dotplot_enrichment(expom, geneset = "GO_BP")
#' }
#'
#' @export
plot_dotplot_enrichment <- function(
    expomicset,
    geneset,
    top_n=5,
    n_per_group=5,
    go_groups=NULL
){
  require(ggplot2)

  if(!"functional_enrichment" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_functional_enrichment` first.")
  }

  go_group_df <- MultiAssayExperiment::metadata(expomicset)$functional_enrichment[[geneset]]

  if(!is.null(go_groups)){
    go_group_df <- go_group_df |>
      dplyr::filter(go_group %in% go_groups)
  }else{
    go_group_df <- go_group_df |>
      dplyr::filter(go_group %in% c(
        go_group_df |>
          dplyr::group_by(go_group) |>
          dplyr::reframe(score=mean(-log10(p.adjust) * Count)) |>
          dplyr::arrange(desc(score)) |>
          dplyr::slice_head(n=top_n) |>
          dplyr::pull(go_group)
      ))
  }

  go_group_df |>
    dplyr::group_by(category,exp_name,go_group) |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n=n_per_group) |>
    dplyr::ungroup() |>
    ggplot(aes(
      x = forcats::fct_reorder(exp_name, -Count),
      y = forcats::fct_reorder(Description, -log10(p.adjust)),
      size = Count,
      color = -log10(p.adjust)
    )) +
    geom_point() +
    ggpubr::theme_pubr(legend = "right",
               base_size = 10) +
    facet_grid(go_group~category, space ="free",scales = "free") +
    ggpubr::rotate_x_text(angle=45)+
    scale_color_gradient(low="thistle1",high="#8E0152") +
    theme(strip.text.x = element_text(face = "bold.italic",angle=90),
          strip.text.y = element_text(face="bold.italic",angle=0))+
    labs(
      x="",
      y="",
      title="",
      size="Number of Genes",
      color=expression(-Log[10](P)),
    )
}
