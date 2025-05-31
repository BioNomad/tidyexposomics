#' Create and store a network from correlation results in a MultiAssayExperiment
#'
#' Constructs an undirected network based on correlation results between exposures
#' and omics data, and optionally stores it in the metadata of a `MultiAssayExperiment` object.
#' Nodes are derived from both variables in the correlation results, and grouped by data type
#' (e.g., "Gene Expression", "Metabolomics", etc.).
#'
#' @param expomicset A `MultiAssayExperiment` object that contains correlation results
#'   in its metadata.
#' @param cor_results Character string indicating which correlation result to use.
#'   Options are `"deg_exp_cor"` (default), `"factor_exp_cor"`, or `"omic_exp_cor"`.
#' @param action Character string indicating whether to `"add"` the resulting graph to
#'   the object's metadata or to `"get"` it directly. Default is `"add"`.
#'
#' @return If `action = "add"`, returns the modified `MultiAssayExperiment` with the
#'   network added to its metadata. If `action = "get"`, returns a list with the graph
#'   object and a summary of the graph structure.
#'
#' @details
#' The function identifies the appropriate correlation results table based on the
#' `cor_results` argument. It then constructs an `igraph` object from the edges
#' (correlations) and derives node metadata from variable types or categories.
#'
#' The resulting graph is undirected and includes group labels for visualization.
#'
#' @import igraph tidygraph MultiAssayExperiment dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a network and add it to the MultiAssayExperiment metadata
#' mae <- run_create_network(mae, cor_results = "omic_exp_cor", action = "add")
#'
#' # Retrieve the network object without adding it
#' net_obj <- run_create_network(mae, cor_results = "deg_exp_cor", action = "get")
#' igraph::plot(net_obj$graph)
#' }

run_create_network <- function(expomicset,
                           cor_results = "deg_exp_cor",
                           action = "add") {

  # Load required libraries
  require(igraph)
  require(tidygraph)

  # Check if data is available
  if(!any(grepl("omics_exposure_deg_correlation|omics_exposure_factor_correlation|omics_exposure_correlation",(MultiAssayExperiment::metadata(expomicset) |> names())))
     ) {
    stop("Please run either `correlate_exposures_degs()` or `correlate_exposures_factors()` or `correlate_exposures_omics()` first.")
  }

  # Switch correlation results based on the input
  cor_res <- switch(cor_results,
                       "deg_exp_cor" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_deg_correlation"]],
                       "factor_exp_cor" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_factor_correlation"]],
                       "omic_exp_cor" = MultiAssayExperiment::metadata(expomicset)[["omics_exposure_correlation"]])

  message("Creating network from correlation results...")

  # Create a graph object from the data frame
  # g <- graph_from_data_frame(cor_res, directed = FALSE)

  # If the correlation results are between omics and exposures:
  if("exp_name" %in% colnames(cor_res)){
    node_tbl <- bind_rows(
      cor_res |> select(name = feature, group = exp_name),
      cor_res |> select(name = exposure, group = category)
    ) |>
      distinct()
  } else if ("var1_type" %in% colnames(cor_res)){
    node_tbl <- bind_rows(
      cor_res |> select(name = var1, group = var1_type),
      cor_res |> select(name = var2, group = var2_type)
    ) |>
      distinct()
  } else {
    stop("Invalid correlation results format.")
  }

  # Create graph with node metadata
  g <- graph_from_data_frame(cor_res, directed = FALSE, vertices = node_tbl)


  # Switch correlation results based on the input
  net_name <- switch(cor_results,
                    "deg_exp_cor" = "omics_exposure_deg_network",
                    "factor_exp_cor" = "omics_exposure_factor_network",
                    "omic_exp_cor" = "omics_exposure_network")

  message("Network Created.")

  if(action=="add"){
    # Save results in metadata
    MultiAssayExperiment::metadata(expomicset)[[net_name]] <- list(
      graph = g,
      summary = .summarize_graph(g)

    )

    # Add analysis steps taken to metadata
    MultiAssayExperiment::metadata(expomicset)$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$steps,
      "run_create_network"
    )

    return(expomicset)
  }else if (action=="get"){
    return(list(
      graph = g,
      summary = .summarize_graph(g)

    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
