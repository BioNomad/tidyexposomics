#' Dotplot of Functional Enrichment Results
#'
#' Generates a dotplot of enriched gene ontology (GO) terms across experiments and categories,
#' optionally annotated with the top genes contributing to each GO group.
#'
#' @param expomicset A `MultiAssayExperiment` object containing functional enrichment results.
#' @param geneset A character string specifying the gene set (e.g., `"GO_BP"` or `"KEGG"`) to extract enrichment results from.
#' @param top_n An integer specifying the number of top GO groups to select (based on enrichment score). Default is `5`.
#' @param n_per_group Number of top GO terms to show within each GO group. Default is `5`.
#' @param add_top_genes Logical; if `TRUE`, appends top genes contributing to each GO group as annotation in facet labels. Default is `TRUE`.
#' @param top_n_genes Integer; number of top genes to extract and display per GO group (used only if `add_top_genes = TRUE`). Default is `5`.
#' @param go_groups Optional character vector of GO group names to include. If `NULL`, the top `top_n` groups are selected automatically.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts enrichment results from `metadata(expomicset)$functional_enrichment[[geneset]]`.
#'   \item Selects top GO groups based on the product of `-log10(p.adjust) * Count`.
#'   \item Displays a dotplot with:
#'     \itemize{
#'       \item x-axis: experimental categories.
#'       \item y-axis: GO term descriptions.
#'       \item dot size: number of genes (`Count`).
#'       \item dot color: significance level (`-log10(p.adjust)`).
#'     }
#'   \item Optionally appends top genes as annotations in facet strip labels.
#' }
#'
#' @return A `ggplot` object displaying the dotplot of enriched GO terms, faceted by GO group and category.
#'
#' @examples
#' \dontrun{
#' plot_dotplot_enrichment(expom, geneset = "GO_BP")
#' plot_dotplot_enrichment(expom, geneset = "KEGG", top_n = 10, add_top_genes = TRUE)
#' }
#'
#' @export

plot_dotplot_enrichment <- function(
    expomicset,
    geneset,
    top_n=5,
    n_per_group=5,
    add_top_genes=TRUE,
    top_n_genes=5,
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

  # go_group_genes_df <- go_group_df |>
  #   (\(df) split(df,df$go_group))() |>
  #   map(~.x |> pull(geneID) |>
  #         (\(x) strsplit(x,"/"))() |>
  #         unlist() |>
  #         table() |>
  #         sort() |>
  #         tail(n=5) |>
  #         names() |>
  #         paste(collapse=",")) |>
  #   as.data.frame() |>
  #   t() |>
  #   as.data.frame() |>
  #   setNames("gene_col") |>
  #   rownames_to_column("go_group")

  go_group_genes_df <- go_group_df |>
    (\(df) split(df, df$go_group))() |>
    purrr::map(~ .x |>
                 pull(geneID) |>
                 (\(x) strsplit(x, "/"))() |>
                 unlist() |>
                 table() |>
                 sort() |>
                 tail(n = top_n_genes) |>
                 names() |>
                 (\(genes) {
                   # Split into groups of 5 and add line breaks
                   gene_chunks <- split(genes, ceiling(seq_along(genes) / 5))
                   paste(sapply(gene_chunks, function(chunk) paste(chunk, collapse = ", ")), collapse = "\n")
                 })()
    ) |>
    as.data.frame() |>
    t() |>
    as.data.frame() |>
    setNames("gene_col") |>
    tibble::rownames_to_column("go_group")

  if(add_top_genes){
    go_group_df <- go_group_df |>
      inner_join(
        go_group_genes_df,
        by = "go_group"
      ) |>
      mutate(go_group = gsub("_"," ", go_group))
  } else{
    go_group_df <- go_group_df |>
      mutate(go_group = gsub("_"," ", go_group))
  }

  go_group_df |>
    mutate(go_group = paste(go_group,gene_col,sep="\n")) |>
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
