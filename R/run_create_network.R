#' Create and store a network from correlation results in a MultiAssayExperiment
#'
#' Constructs an undirected network based on correlation results between exposures
#' and features (e.g., DEGs, omics, latent factors, exposures) and optionally stores it
#' in the metadata of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param feature_type One of "degs", "omics", "factors", or "exposures". Used to select correlation results.
#' @param action "add" to store the graph in metadata, or "get" to return it.
#'
#' @return If `action = "add"`, returns modified `MultiAssayExperiment`. Otherwise, returns a list with `graph` and `summary`.
#' @export
run_create_network <- function(expomicset,
                               feature_type = c("degs", "omics", "factors", "exposures"),
                               action = c("add", "get")) {
  require(igraph)
  require(tidygraph)

  feature_type <- match.arg(feature_type)
  action <- match.arg(action)

  cor_data <- MultiAssayExperiment::metadata(expomicset)$correlation[[feature_type]]

  if (is.null(cor_data)) {
    stop(paste0("No correlation data found for feature_type = '", feature_type, "'"))
  }

  message("Creating network from correlation results...")

  # Detect variable structure
  if (all(c("feature", "exposure", "exp_name", "category") %in% colnames(cor_data))) {
    node_tbl <- dplyr::bind_rows(
      dplyr::select(cor_data, name = feature, group = exp_name),
      dplyr::select(cor_data, name = exposure, group = category)
    ) |> dplyr::distinct()

    edge_df <- dplyr::select(cor_data, from = exposure, to = feature, correlation, FDR)

  } else if (all(c("var1", "var2", "var1_type", "var2_type") %in% colnames(cor_data))) {
    node_tbl <- dplyr::bind_rows(
      dplyr::select(cor_data, name = var1, group = var1_type),
      dplyr::select(cor_data, name = var2, group = var2_type)
    ) |> dplyr::distinct()

    edge_df <- dplyr::select(cor_data, from = var1, to = var2, correlation, FDR)

  } else {
    stop("Unrecognized correlation format.")
  }

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_tbl)
  net_summary <- .summarize_graph(g)

  net_name <- paste0("network_", feature_type)

  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)[["network"]][[net_name]] <- list(
      graph = g,
      summary = net_summary
    )

    # Add analysis step to metadata
    step_record <- list(
      run_create_network = list(
        timestamp = Sys.time(),
        params = list(
          feature_type = feature_type
        ),
        notes = paste0(
          "Created undirected network from correlation results for feature_type = '",
          feature_type,
          "'.")
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    message("Network added to metadata as: ", net_name)
    return(expomicset)
  } else {
    return(list(graph = g, summary = net_summary))
  }
}
