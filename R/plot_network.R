#' Plot Network Graph of Omics-Exposure Associations
#'
#' Visualizes a tidygraph network stored in the `MultiAssayExperiment` metadata using `ggraph`.
#'
#' @param expomicset A `MultiAssayExperiment` object containing network metadata (e.g., `"omics_exposure_network"`).
#' @param network A character string specifying which network to plot. Options are:
#' \describe{
#'   \item{"omics_exposure_network"}{The full correlation network.}
#'   \item{"omics_exposure_deg_network"}{Subset of the network for DEGs.}
#'   \item{"omics_exposure_factor_network"}{Subset of the network involving factor features.}
#' }
#' @param include_stats Logical; whether to compute and visualize centrality statistics. Default is `TRUE`.
#' @param nodes_to_include A character vector of node names to retain. If `NULL`, include all nodes.
#' @param centrality_thresh Numeric threshold to retain only nodes above a minimum centrality value.
#' @param top_n_nodes Integer; retain only the top N most central nodes.
#' @param cor_thresh Numeric; threshold for filtering edges by absolute correlation.
#' @param label Logical; whether to display node labels. Default is `FALSE`.
#' @param label_top_n Integer; number of top central nodes to label if `label = TRUE` and `nodes_to_label` is `NULL`. Default is `5`.
#' @param nodes_to_label A character vector of node names to label. Overrides `label_top_n` if provided.
#' @param facet_var Optional column name (in the node data) to use for faceting the layout.
#' @param foreground Color of the label highlight. Default is `"steelblue"`.
#' @param fg_text_colour Text color for node labels. Default is `"grey25"`.
#' @param node_colors Optional named vector for manually setting node colors by group.
#' @param node_color_var A column name in node metadata used for color mapping (e.g., `"type"` or `"category"`).
#' @param alpha Transparency for nodes. Default is `0.5`.
#' @param size_lab Label for the node size legend. Default is `"Centrality"`.
#' @param color_lab Label for the node color legend. Default is `"Group"`.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Selects a stored network from `metadata(expomicset)`.
#'   \item Applies optional node and edge filters (e.g., correlation threshold, centrality, node list).
#'   \item Prunes unconnected nodes (nodes not involved in any remaining edges).
#'   \item Computes node centrality for sizing or filtering if requested.
#'   \item Generates a `ggraph` layout using `.build_ggraph_plot()`.
#' }
#' Node color and label aesthetics are customizable. Labeling can be automatic (e.g., top 5 by centrality) or manual via `nodes_to_label`.
#'
#' @return A `ggraph` object for plotting.
#'
#' @examples
#' \dontrun{
#' plot_network(
#'   expomicset,
#'   network = "omics_exposure_network",
#'   cor_thresh = 0.4,
#'   top_n_nodes = 100,
#'   label = TRUE
#' )
#' }
#'
#' @export

plot_network <- function(expomicset,
                         network,
                         include_stats = TRUE,
                         nodes_to_include = NULL,
                         centrality_thresh = NULL,
                         top_n_nodes = NULL,
                         cor_thresh = NULL,
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

  if(!is.null(cor_thresh)){
    message("Filtering nodes based on correlation threshold...")
    g <- g |>
      activate(edges) |>
      filter(abs(correlation) > cor_thresh)
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

  # Filter the graph before plotting to make sure that
  # only nodes included in edges are present
  # Keep only nodes involved in at least one edge
  used_node_indices <- g |>
    activate(edges) |>
    as_tibble() |>
    select(from, to) |>
    unlist() |>
    unique()

  g <- g |>
    activate(nodes) |>
    mutate(node_index = row_number()) |>
    filter(node_index %in% used_node_indices) |>
    select(-node_index)

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
