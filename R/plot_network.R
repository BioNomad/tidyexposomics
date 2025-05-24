#' Plot a tidygraph network from a MultiAssayExperiment object
#'
#' This function visualizes a network stored in the metadata of a `MultiAssayExperiment`
#' object using `ggraph`. It supports node filtering, labeling, faceting, and visual annotation
#' based on node attributes such as centrality or metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object that contains a network in its metadata.
#' @param network A character string indicating which network to extract. Options:
#'   `"omics_exposure_deg_network"`, `"omics_exposure_factor_network"`, or `"omics_exposure_network"`.
#' @param include_stats Logical; if `TRUE`, adds a caption summarizing graph statistics.
#' @param nodes_to_include Optional character vector of node names to retain in the plot.
#' @param centrality_thresh Optional numeric threshold to filter nodes by degree centrality.
#' @param top_n_nodes Optional integer; keeps the top N nodes by centrality.
#' @param label Logical; if `TRUE`, adds node labels.
#' @param label_top_n Integer; number of top central nodes to label if `nodes_to_label` is not provided.
#' @param nodes_to_label Optional character vector; names of nodes to explicitly label.
#' @param facet_var Optional name of a node attribute to facet the plot by (e.g., `group`, `type`).
#' @param foreground Background fill color for facet strip themes.
#' @param fg_text_colour Text color for graph theme.
#' @param node_colors Optional named vector of colors to apply to node groups (used with `node_color_var`).
#' @param node_color_var Optional name of a node attribute to color nodes by (e.g., `type`, `module`).
#' @param alpha Numeric; transparency level for facet strip backgrounds.
#' @param size_lab Label for the node size legend (default: `"Centrality"`).
#' @param color_lab Label for the node color legend (default: `"Group"`).
#'
#' @return A `ggplot` object representing the plotted network.
#' @import ggraph tidygraph igraph
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume `mae` is a MultiAssayExperiment with a network
#' plot_network(mae, network = "omics_exposure_network")
#' }

plot_network <- function(expomicset,
                         network,
                         include_stats = TRUE,
                         nodes_to_include = NULL,
                         centrality_thresh = NULL,
                         top_n_nodes = NULL,
                         label = FALSE,
                         label_top_n = 5,
                         nodes_to_label = NULL,
                         facet_var = NULL,
                         foreground = 'steelblue',
                         fg_text_colour = 'grey25',
                         node_colors = NULL,
                         node_color_var = NULL,
                         alpha = 0.5,
                         size_lab = "Centrality",
                         color_lab = "Group"
                         ){
  # Load required libraries
  require(ggraph)
  require(tidygraph)

  # Check to see that networks are present in the mae
  if(!any(grepl("network",(MultiAssayExperiment::metadata(expomicset) |> names())))
  ) {
    stop("Please run either `run_create_network()` first.")
  }

  # Switch the graph based on which network user choses
  g <- switch(network,
                     "omics_exposure_deg_network" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_deg_network"]][["graph"]],
                     "omics_exposure_factor_network" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_factor_network"]][["graph"]],
                     "omics_exposure_network" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_network"]][["graph"]])

  message("Extracting graph...")

  # Ensure graph is a tidygraph object
  g <- g |>
    tidygraph::as_tbl_graph()

  # Filter based on whether or not the user wants certain nodes
  if(!is.null(nodes_to_include)){
    message("Filtering nodes...")
    g <- g |>
      activate(nodes) |>
      filter(name %in% nodes_to_include)
  }

  # Filter based on centrality threshold
  if(!is.null(centrality_thresh)){
    message("Filtering nodes based on centrality threshold...")
    g <- g |>
      activate(nodes) |>
      mutate(centrality=tidygraph::centrality_degree()) |>
      filter(centrality>centrality_thresh)
  }

  # Filter top n nodes based on centrality
  if(!is.null(top_n_nodes)){
    message("Filtering top ",top_n_nodes," nodes based on centrality...")
    g <- g |>
      activate(nodes) |>
      mutate(centrality=tidygraph::centrality_degree()) |>
      arrange(desc(centrality)) |>
      slice_head(n=top_n_nodes)
  }

  # Create a column of nodes to label based on user input
  if(!is.null(nodes_to_label) & isTRUE(label)){
    g <- g |>
      activate(nodes) |>
      mutate(label=ifelse(name %in% nodes_to_label, name, NA))

  } else if(isTRUE(label)){
    # If the user wants to label nodes, but does not provide input
    # label the top 5 nodes based on centrality
    g <- g |>
      activate(nodes) |>
      mutate(centrality=tidygraph::centrality_degree()) |>
      arrange(desc(centrality)) |>
      mutate(label = ifelse(dplyr::row_number() <= label_top_n, name, NA))
  }

  # User confirmation for large networks
  if (igraph::gorder(g) > 500) {
    user_input <- readline("The network has more than 500 nodes. Plot anyway? (y/n): ")
    if (tolower(user_input) != "y") stop("Exiting.")
  }

  # Build and return plot using helper
  g_plot <- .build_ggraph_plot(
    g = g,
    node_color_var = node_color_var,
    label = label,
    facet_var = facet_var,
    include_stats = include_stats,
    fg_text_colour = fg_text_colour,
    foreground = foreground,
    alpha = alpha,
    color_lab = color_lab,
    node_colors = node_colors,
    size_lab = size_lab
  )

  return(g_plot)
}
