#' Plot Network Graph of Features or Exposures
#'
#' Visualizes network structures created by `run_create_network()` from
#' the metadata of a `MultiAssayExperiment` object.
#' Nodes can represent features (e.g., genes or factors) or exposures,
#' and edges represent correlations or shared connections.
#'
#' @param expomicset A `MultiAssayExperiment` object containing network results
#' in metadata.
#' @param network Character string specifying the network type.
#' One of `"degs"`, `"omics"`, `"factors"`,
#'   `"factor_features"`, `"exposures"`, `"degs_feature_cor"`,
#'   `"omics_feature_cor"`, `"factor_features_feature_cor"`.
#' @param include_stats Logical; if `TRUE`, include edge weights and
#' node centrality metrics in the plot aesthetics. Default is `TRUE`.
#' @param nodes_to_include Optional character vector of node names to
#' include (subset of `name`).
#' @param centrality_thresh Optional numeric threshold to filter nodes
#' by centrality degree.
#' @param top_n_nodes Optional integer to keep only the top N nodes by
#' centrality.
#' @param cor_thresh Optional numeric threshold to filter edges by
#' minimum absolute correlation.
#' @param label Logical; whether to label nodes.
#' If `TRUE`, top nodes will be labeled.
#' @param label_top_n Integer; number of top-centrality nodes to
#' label if `label = TRUE`. Default is `5`.
#' @param nodes_to_label Optional character vector of specific nodes to label.
#' @param facet_var Optional node metadata column to facet the
#' network layout by.
#' @param foreground Color for node outlines and edge borders.
#' Default is `"steelblue"`.
#' @param fg_text_colour Color of node label text. Default is `"grey25"`.
#' @param node_colors Optional named vector of colors for node groups.
#' @param node_color_var Optional node attribute used for node color mapping.
#' @param alpha Alpha transparency for nodes and edges. Default is `0.5`.
#' @param size_lab Legend title for node size (typically centrality).
#'  Default is `"Centrality"`.
#' @param color_lab Legend title for node color group. Default is `"Group"`.
#'
#' @return A `ggraph` plot object.
#'
#' @details
#' This function retrieves the stored graph object and optionally filters or
#'  labels nodes based on:
#' centrality, correlation, user input, or group-specific attributes.
#'  It supports layout faceting,
#' custom color mappings, and highlights highly central nodes.
#'
#' Large graphs (> 500 nodes) will prompt the user before plotting.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # perform correlation analyses
#' # correlate with exposures
#' mae <- mae |>
#'     run_correlation(
#'         feature_type = "omics",
#'         variable_map = mae |>
#'             pivot_feature() |>
#'             dplyr::select(
#'                 variable = .feature,
#'                 exp_name = .exp_name
#'             ),
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' # create the networks
#' mae <- mae |>
#'     run_create_network(
#'         feature_type = "omics",
#'         action = "add"
#'     )
#'
#' # plot the network
#' network_p <- mae |>
#'     plot_network(
#'         network = "omics"
#'     )
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom tidygraph as_tbl_graph activate filter mutate centrality_degree
#'   arrange slice_head select as_tibble
#' @importFrom igraph gorder
#' @export
plot_network <- function(
    expomicset,
    network = c(
        "degs",
        "omics",
        "factors",
        "factor_features",
        "exposures",
        "degs_feature_cor",
        "omics_feature_cor",
        "factor_features_feature_cor"
    ),
    include_stats = TRUE,
    nodes_to_include = NULL,
    centrality_thresh = NULL,
    top_n_nodes = NULL,
    cor_thresh = NULL,
    label = FALSE,
    label_top_n = 5,
    nodes_to_label = NULL,
    facet_var = NULL,
    foreground = "steelblue",
    fg_text_colour = "grey25",
    node_colors = NULL,
    node_color_var = NULL,
    alpha = 0.5,
    size_lab = "Centrality",
    color_lab = "Group") {
    # require(ggraph)
    # require(tidygraph)

    network <- match.arg(network)

    net_key <- paste0("network_", network)
    net_obj <- MultiAssayExperiment::metadata(expomicset)$network[[net_key]]

    if (is.null(net_obj)) {
        stop(
            "No network found for `feature_type = '",
            network,
            "'. Please run `run_create_network()` first."
        )
    }

    message("Extracting graph.")
    g <- tidygraph::as_tbl_graph(net_obj$graph)

    if (!is.null(nodes_to_include)) {
        g <- g |>
            activate(nodes) |>
            filter(name %in% nodes_to_include)
    }

    if (!is.null(cor_thresh)) {
        g <- g |>
            activate(edges) |>
            filter(abs(correlation) > cor_thresh)
    }

    if (!is.null(centrality_thresh)) {
        g <- g |>
            activate(nodes) |>
            mutate(centrality = centrality_degree()) |>
            filter(centrality > centrality_thresh)
    }

    if (!is.null(top_n_nodes)) {
        g <- g |>
            activate(nodes) |>
            mutate(centrality = centrality_degree()) |>
            arrange(desc(centrality)) |>
            slice_head(n = top_n_nodes)
    }

    if (!is.null(nodes_to_label) && isTRUE(label)) {
        g <- g |>
            activate(nodes) |>
            mutate(label = ifelse(name %in% nodes_to_label, name, NA))
    } else if (isTRUE(label)) {
        g <- g |>
            activate(nodes) |>
            mutate(centrality = centrality_degree()) |>
            arrange(desc(centrality)) |>
            mutate(label = ifelse(row_number() <= label_top_n, name, NA))
    }

    # Prune nodes without edges
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

    # Warn if plotting very large networks
    if (igraph::gorder(g) > 500) {
        user_input <- readline("The network has more than 500 nodes.
        Plot anyway? (y/n): ")
        if (tolower(user_input) != "y") stop("Exiting.")
    }

    # Build plot using internal helper
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
