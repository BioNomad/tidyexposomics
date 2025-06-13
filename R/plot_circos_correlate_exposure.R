#' Plot Circular Network of Exposure Correlations
#'
#' Generates a circular network plot to visualize significant exposure-exposure correlations from a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure correlation results in `metadata(expomicset)$exposure_correlation$correlation_table`.
#' @param exposure_cols Optional character vector specifying which exposure variables to include. If `NULL`, all are used.
#' @param corr_threshold Optional numeric value specifying the minimum absolute correlation for inclusion. If `NULL`, no thresholding is applied.
#' @param annotation_colors Optional named vector of colors for exposure categories. Names must match exposure category labels.
#' @param low Color used for low correlations in the edge color gradient. Default is "#006666".
#' @param mid Color used for zero correlation in the edge color gradient. Default is "white".
#' @param high Color used for high correlations in the edge color gradient. Default is "#8E0152".
#' @param midpoint Numeric value specifying the midpoint of the color gradient. Default is 0.
#'
#' @details
#' This function builds a correlation network of exposures, where edges represent significant correlations between exposures (filtered by `corr_threshold`)
#' and nodes represent individual exposure variables. Nodes are colored by category (from `metadata(expomicset)$var_info`) and arranged in a circular layout.
#'
#' Correlation results must have been previously calculated and stored in `metadata(expomicset)$exposure_correlation$correlation_table`
#' using a function like `correlate_exposures()`. Exposure metadata is retrieved from `metadata(expomicset)$var_info`.
#'
#' @return A `ggplot` object created with `ggraph`, visualizing the exposure correlation network in a circular layout.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_circos_correlate_exposure(expom)
#'
#' # With custom threshold and subset of exposures
#' plot_circos_correlate_exposure(
#'   expom,
#'   exposure_cols = c("pm25", "no2", "o3"),
#'   corr_threshold = 0.4
#' )
#' }
#' @export
plot_circos_correlate_exposure <- function(
    expomicset,
    exposure_cols = NULL,
    corr_threshold = NULL,
    annotation_colors = NULL,
    low = "#006666",
    mid = "white",
    high = "#8E0152",
    midpoint = 0) {

  require(ggplot2)
  require(ggraph)

  # Check that exposure correlation has been run
  if(!("exposure_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))) ) {
    stop("Exposure correlation has not been run. Please run `correlate_exposures()` first.")
  }

  correlation_df <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("exposure_correlation") |>
    purrr::pluck("correlation_table")

  if(!is.null(exposure_cols)){
    # Filter correlation_df based on exposure_cols
    correlation_df <- correlation_df |>
      dplyr::filter(var1 %in% exposure_cols &
                      var2 %in% exposure_cols)
  }

  if(!is.null(corr_threshold)) {
    # Filter based on correlation threshold
    correlation_df <- correlation_df |>
      dplyr::filter(abs(correlation) >= corr_threshold)

    # Remove unnecessary rows if no correlations remain
    remaining_vars <- unique(c(correlation_df$var1, correlation_df$var2))

    correlation_df <- correlation_df |>
      dplyr::filter(var1 %in% remaining_vars & var2 %in% remaining_vars)
  }

  # Unified color palette
  all_categories <- unique(c(correlation_df$category_1,
                             correlation_df$category_2)) |>
    as.character() |>
    unique() |>
    na.omit()

  # Check if annotation_colors is provided
  if(!is.null(annotation_colors)){
    cat_colors <- annotation_colors

    # Check if the length of annotation_colors matches the number of unique categories
    if(length(cat_colors) != length(all_categories)) {
      stop("Length of annotation_colors must match number of unique categories in the data.")
    }

  } else{
    cat_colors <- setNames(tidy_exp_pal[1:length(all_categories)], all_categories)
  }

  # Vertex Information
  vertex_df <- MultiAssayExperiment::metadata(expomicset)$var_info |>
    dplyr::filter(variable %in% c(correlation_df$var1,
                                  correlation_df$var2))

  # Create igraph object
  graph <- igraph::graph_from_data_frame(correlation_df,
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
      low=low,
      mid=mid,
      high=high,
      midpoint=midpoint)+
    ggplot2::coord_fixed(xlim = c(-2, 2),
                ylim = c(-2, 2)) +
    ggplot2::guides(edge_width = "none")+
    ggplot2::labs(
      edge_color = "Correlation",
      fill = "Category"
    )
}



