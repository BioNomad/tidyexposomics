#' Plot Circular Network of Exposures Sharing Common Features
#'
#' Generates a circular correlation network plot to visualize exposures that share a high number of associated features.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use.
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#' @param cutoff An integer specifying the minimum number of shared features required to display an exposure-exposure connection. Default is `10`.
#'
#' @details
#' This function identifies exposures that share a significant number of correlated features.
#' It constructs a circular network visualization using `ggraph`, where:
#' - Nodes represent exposures, with fill color indicating their category.
#' - Edges represent the number of shared features, with thicker and more intense edges for stronger overlaps.
#'
#' The input correlation data is extracted from `metadata(expomicset)`, and only exposure pairs sharing more than `cutoff` features are included.
#'
#' @return A `ggraph` plot displaying a circular network of exposures sharing common features.
#'
#' @examples
#' \dontrun{
#' plot_circos_exposure_shared_features(expom)
#' }
#'
#' @export
plot_circos_exposure_shared_features <- function(
    expomicset,
    geneset = "degs",
    cutoff = 10) {

  require(ggplot2)
  require(ggraph)


  if(geneset=="degs"){
    # Check if the metadata contains the required data
    if(!"omics_exposure_deg_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation

  }else if(geneset=="factors"){
    # Check if the metadata contains the required data
    if(!"omics_exposure_factor_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }

  # Extract correlated features
  df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation

  # Filter the data to only include features with correlation above the cutoff
  exp_feature_cor_lst <- split(
    exp_feature_cor_df,
    exp_feature_cor_df$exposure) |>
    purrr::map(~ .x |>
          dplyr::pull(feature) |>
          unique())

  # Get pairwise overlaps
  exp_shared_feature_df <- .get_pairwise_overlaps(exp_feature_cor_lst) |>
    dplyr::filter(num_shared > cutoff) |>
    dplyr::arrange(num_shared)

  # Vertex Information
  vertex_df <- MultiAssayExperiment::metadata(expomicset)$var_info %>%
    dplyr::filter(variable %in% c(
      exp_shared_feature_df$source,
      exp_shared_feature_df$target))

  # Create igraph object
  graph <- igraph::graph_from_data_frame(exp_shared_feature_df,
                                 vertices = vertex_df,
                                 directed = FALSE)

  # Plot using ggraph
  ggraph::ggraph(graph,
         layout = "linear",
         circular = TRUE) +
    ggraph::geom_edge_arc(aes(
      color = num_shared,
      width = num_shared)) +
    ggraph::geom_node_text(aes(
      label = name,
      fontface = "bold.italic",
      angle = ggraph::node_angle(x, y)),
      hjust = "outward",
      check_overlap = TRUE) +
    ggraph::geom_node_point(shape = 21,
                    size = 4,
                    alpha=0.8,
                    aes(fill = category)) +
    ggraph::theme_graph() +
    ggsci::scale_fill_aaas()+
    ggraph::scale_edge_color_gradient2(
      low="#006666",
      mid="white",
      high="#8E0152",
      midpoint = mean(exp_shared_feature_df$num_shared) )+
    coord_fixed(xlim = c(-2, 2),
                ylim = c(-2, 2)) +
    #theme(text = element_text(face = "bold"))+
    guides(edge_width = "none")+
    labs(
      edge_width="Number of Shared Features",
      edge_color="Number of Shared Features",
      fill="Category"
    )
}
