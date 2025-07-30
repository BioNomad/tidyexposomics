#' Create Correlation Network from Feature Data
#'
#' Constructs an undirected feature–feature or feature–exposure
#' correlation network from correlation results stored in a
#' `MultiAssayExperiment` object. The function
#' supports multiple correlation formats depending on `feature_type`,
#' and stores or returns an `igraph` object with associated node
#' and edge metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object containing correlation
#' results in metadata.
#' @param feature_type Type of correlation result to convert to a network.
#'  One of:
#'   `"degs"`, `"omics"`, `"factors"`, `"factor_features"`, `"exposures"`,
#'   `"degs_feature_cor"`, `"omics_feature_cor"`,
#'    or `"factor_features_feature_cor"`.
#' @param action Whether to `"add"` the network to the object or
#' `"get"` it as a list.
#'
#' @return If `action = "add"`, returns the updated `MultiAssayExperiment`
#' with a new `network` entry in metadata.
#'         If `action = "get"`, returns a list with `graph`
#'         (an `igraph` object) and `summary` (a tibble).
#'
#' @details The function detects the appropriate edge and node structure
#' based on column names in the correlation results.
#' Edge weights are based on correlation coefficients and include FDR values.
#'
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 10,
#'   return_mae=TRUE
#'   )
#'
#' # perform correlation analyses
#' # correlate with exposures
#' mae <- mae |>
#'   run_correlation(
#'   feature_type = "omics",
#'   variable_map = mae |>
#'     pivot_feature() |>
#'     dplyr::select(variable=.feature,
#'     exp_name=.exp_name),
#'   exposure_cols = c("exposure_pm25","exposure_no2","age","bmi")
#'   ) |>
#'   # correlate omics features with themselves
#'   run_correlation(
#'   feature_type = "omics",
#'   variable_map = mae |>
#'     pivot_feature() |>
#'     dplyr::select(variable=.feature,
#'     exp_name=.exp_name),
#'   feature_cors = TRUE,
#'   exposure_cols = c("exposure_pm25","exposure_no2","age","bmi")
#'   )
#'
#' # create the networks
#' mae   <- mae |>
#' run_create_network(
#'   feature_type = "omics_feature_cor",
#'   action = "add") |>
#'   run_create_network(
#'    feature_type = "omics",
#'    action = "add")
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom dplyr bind_rows select distinct
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @export
run_create_network <- function(expomicset,
                               feature_type = c("degs",
                                                "omics",
                                                "factors",
                                                "factor_features",
                                                "exposures",
                                                "degs_feature_cor",
                                                "omics_feature_cor",
                                                "factor_features_feature_cor"),
                               action = c("add",
                                          "get")) {
  # require(igraph)
  # require(tidygraph)

  feature_type <- match.arg(feature_type)
  action <- match.arg(action)

  all_metadata <- MultiAssayExperiment::metadata(expomicset)
  cor_data <- all_metadata$correlation[[feature_type]]


  if (is.null(cor_data)) {
    stop(sprintf(
      "No correlation data found for feature_type = '%s'",
      feature_type))
  }


  message("Creating network from correlation results.")

  # Detect variable structure
  if (all(c("feature",
            "exposure",
            "exp_name",
            "category") %in% colnames(cor_data))) {
    node_tbl <- dplyr::bind_rows(
      dplyr::select(cor_data, name = feature, group = exp_name),
      dplyr::select(cor_data, name = exposure, group = category)
    ) |> dplyr::distinct()

    edge_df <- dplyr::select(cor_data,
                             from = exposure,
                             to = feature,
                             correlation, FDR)

  } else if (all(c("var1",
                   "var2",
                   "var1_type",
                   "var2_type") %in% colnames(cor_data))) {
    node_tbl <- dplyr::bind_rows(
      dplyr::select(cor_data, name = var1, group = var1_type),
      dplyr::select(cor_data, name = var2, group = var2_type)
    ) |> dplyr::distinct()

    edge_df <- dplyr::select(cor_data,
                             from = var1,
                             to = var2,
                             correlation, FDR)

  } else if (all(c("var1",
                   "var2",
                   "exp_name_1",
                   "exp_name_2") %in% colnames(cor_data))){
    node_tbl <- dplyr::bind_rows(
      dplyr::select(cor_data, name = var1, group = exp_name_1),
      dplyr::select(cor_data, name = var2, group = exp_name_2)
    ) |> dplyr::distinct()

    edge_df <- dplyr::select(cor_data,
                             from = var1,
                             to = var2,
                             correlation,
                             FDR)
  } else {
    stop("Unrecognized correlation format.")
  }

  g <- igraph::graph_from_data_frame(edge_df,
                                     directed = FALSE,
                                     vertices = node_tbl)
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
          "Created undirected network from correlation results for'",
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
