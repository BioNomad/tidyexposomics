% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_network.R
\name{plot_network}
\alias{plot_network}
\title{Plot Network Graph of Omics-Exposure Associations}
\usage{
plot_network(
  expomicset,
  network = c("degs", "omics", "factors", "exposures"),
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
  color_lab = "Group"
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing network metadata (e.g., \code{"omics_exposure_network"}).}

\item{network}{A character string specifying which network to plot. Options are:
\describe{
\item{"omics_exposure_network"}{The full correlation network.}
\item{"omics_exposure_deg_network"}{Subset of the network for DEGs.}
\item{"omics_exposure_factor_network"}{Subset of the network involving factor features.}
}}

\item{include_stats}{Logical; whether to compute and visualize centrality statistics. Default is \code{TRUE}.}

\item{nodes_to_include}{A character vector of node names to retain. If \code{NULL}, include all nodes.}

\item{centrality_thresh}{Numeric threshold to retain only nodes above a minimum centrality value.}

\item{top_n_nodes}{Integer; retain only the top N most central nodes.}

\item{cor_thresh}{Numeric; threshold for filtering edges by absolute correlation.}

\item{label}{Logical; whether to display node labels. Default is \code{FALSE}.}

\item{label_top_n}{Integer; number of top central nodes to label if \code{label = TRUE} and \code{nodes_to_label} is \code{NULL}. Default is \code{5}.}

\item{nodes_to_label}{A character vector of node names to label. Overrides \code{label_top_n} if provided.}

\item{facet_var}{Optional column name (in the node data) to use for faceting the layout.}

\item{foreground}{Color of the label highlight. Default is \code{"steelblue"}.}

\item{fg_text_colour}{Text color for node labels. Default is \code{"grey25"}.}

\item{node_colors}{Optional named vector for manually setting node colors by group.}

\item{node_color_var}{A column name in node metadata used for color mapping (e.g., \code{"type"} or \code{"category"}).}

\item{alpha}{Transparency for nodes. Default is \code{0.5}.}

\item{size_lab}{Label for the node size legend. Default is \code{"Centrality"}.}

\item{color_lab}{Label for the node color legend. Default is \code{"Group"}.}
}
\value{
A \code{ggraph} object for plotting.
}
\description{
Visualizes a tidygraph network stored in the \code{MultiAssayExperiment} metadata using \code{ggraph}.
}
\details{
This function:
\itemize{
\item Selects a stored network from \code{metadata(expomicset)}.
\item Applies optional node and edge filters (e.g., correlation threshold, centrality, node list).
\item Prunes unconnected nodes (nodes not involved in any remaining edges).
\item Computes node centrality for sizing or filtering if requested.
\item Generates a \code{ggraph} layout using \code{.build_ggraph_plot()}.
}
Node color and label aesthetics are customizable. Labeling can be automatic (e.g., top 5 by centrality) or manual via \code{nodes_to_label}.
}
\examples{
\dontrun{
plot_network(
  expomicset,
  network = "omics_exposure_network",
  cor_thresh = 0.4,
  top_n_nodes = 100,
  label = TRUE
)
}

}
