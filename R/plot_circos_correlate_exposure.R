#' Plot Circular Network of Exposure Correlations
#'
#' Generates a circular correlation network plot to visualize significant exposure-exposure correlations.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure correlation results.
#' @param correlation_cutoff A numeric value specifying the minimum absolute correlation required for inclusion in the plot. Default is `0.3`.
#'
#' @details
#' This function extracts exposure correlation results from `metadata(expomicset)$exposure_correlation$filtered_table`,
#' filters edges based on `correlation_cutoff`, and constructs a circular network visualization using `ggraph`.
#' Nodes represent exposures, with fill color indicating exposure category, while edges represent significant correlations.
#'
#' @return A `ggraph` plot displaying a circular correlation network of exposures.
#'
#' @examples
#' \dontrun{
#' plot_circos_correlate_exposure(expom)
#' }
#'
#' @export
plot_circos_correlate_exposure <- function(
    expomicset,
    correlation_cutoff = 0.3) {

  require(ggplot2)
  require(ggraph)

  # Extract and filter relevant data
  correlation_data <- MultiAssayExperiment::metadata(expomicset)$exposure_correlation$filtered_table %>%
    dplyr::filter(abs_correlation >= correlation_cutoff)

  # Vertex Information
  vertex_df <- MultiAssayExperiment::metadata(expomicset)$var_info %>%
    dplyr::filter(variable %in% c(correlation_data$var1, correlation_data$var2))

  # Create igraph object
  graph <- igraph::graph_from_data_frame(correlation_data,
                                 vertices = vertex_df,
                                 directed = FALSE)

  # Plot using ggraph
  ggraph::ggraph(graph,
         layout = "linear",
         circular = TRUE) +
    ggraph::geom_edge_arc(aes(
      color = correlation,
      width = abs_correlation)) +
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
    ggsci::scale_fill_npg()+
    ggraph::scale_edge_color_gradient2(
      low="#006666",
      mid="white",
      high="#8E0152",
      midpoint=0)+
    ggplot2::coord_fixed(xlim = c(-2, 2),
                ylim = c(-2, 2)) +
    ggplot2::guides(edge_width = "none")+
    ggplot2::labs(
      edge_color = "Correlation",
      fill = "Category"
    )
}



