#' Calculate Exposure Impact from Feature-Exposure Correlation Networks
#'
#' Generalized centrality-based exposure impact analysis using DEG, omics,
#' or factor features.
#'
#' @param expomicset A `MultiAssayExperiment` object with correlation
#' and network metadata.
#' @param feature_type One of `"degs"`, `"omics"`, or `"factor_features"`.
#' @param pval_col Column in differential abundance results to filter DEGs.
#' Default = `"adj.P.Val"`.
#' @param pval_thresh DEG p-value threshold. Ignored unless
#' `feature_type == "degs"`.
#' @param action Either `"add"` (store in metadata) or `"get"` (return list).
#'
#' @return Either an updated MultiAssayExperiment (if action = "add") or a list.
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
#'     ) |>
#'     run_correlation(
#'         feature_type = "omics",
#'         variable_map = mae |>
#'             pivot_feature() |>
#'             dplyr::select(
#'                 variable = .feature,
#'                 exp_name = .exp_name
#'             ),
#'         feature_cors = TRUE,
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' # create the networks
#' mae <- mae |>
#'     run_create_network(
#'         feature_type = "omics_feature_cor",
#'         action = "add"
#'     ) |>
#'     run_create_network(
#'         feature_type = "omics",
#'         action = "add"
#'     )
#'
#' # perform impact analysis
#' mae <- mae |>
#'     run_exposure_impact(
#'         feature_type = "omics"
#'     )
#'
#' @export
run_exposure_impact <- function(
    expomicset,
    feature_type = c(
        "degs",
        "omics",
        "factor_features"
    ),
    pval_col = "adj.P.Val",
    pval_thresh = 0.1,
    action = c("add", "get")) {
    .check_suggested(pkg = "tidygraph")
    feature_type <- match.arg(feature_type)
    action <- match.arg(action)
    md <- MultiAssayExperiment::metadata(expomicset)

    # Validate required correlation and network
    if (is.null(md$correlation[[feature_type]])) {
        stop(
            "Correlation results missing for feature_type = ",
            feature_type,
            ". Please run run_correlation(..., feature_type = '",
            feature_type, "')."
        )
    }
    if (is.null(md$network[[paste0("network_", feature_type)]])) {
        stop(
            "Network missing for feature_type = ", feature_type,
            ". Please run run_create_network(..., feature_type = '",
            feature_type, "')."
        )
    }

    # Get correlation and network
    correlation_df <- md$correlation[[feature_type]]
    graph_list_name <- paste0("network_", feature_type, "_feature_cor")
    net_graph <- md$network[[graph_list_name]]$graph

    # Compute centrality for all nodes
    # node_centrality <- net_graph |>
    #   tidygraph::as_tbl_graph() |>
    #   tidygraph::activate(nodes) |>
    #   dplyr::mutate(centrality = tidygraph::centrality_degree()) |>
    #   dplyr::as_tibble()

    node_centrality <- net_graph |>
        tidygraph::as_tbl_graph() |>
        tidygraph::activate(nodes) |>
        dplyr::mutate(
            centrality_degree = tidygraph::centrality_degree(),
            centrality_eigen = tidygraph::centrality_eigen(),
            centrality_closeness = tidygraph::centrality_closeness(),
            centrality_betweenness = tidygraph::centrality_betweenness()
        ) |>
        dplyr::as_tibble()

    # Merge with correlation edges (exposure-feature)
    if (feature_type == "degs") {
        correlation_df <- correlation_df |>
            dplyr::rename(
                feature = feature,
                exposure = exposure
            )
    } else if ("var1" %in% names(correlation_df)) {
        # Handle generalized format
        correlation_df <- correlation_df |>
            dplyr::filter(
                var1_type == "exposure",
                var2_type != "exposure"
            ) |>
            dplyr::rename(
                exposure = var1,
                feature = var2
            )
    }

    impact_df <- correlation_df |>
        dplyr::left_join(node_centrality, by = c("feature" = "name")) |>
        dplyr::filter(!is.na(centrality_degree)) |>
        dplyr::group_by(exposure) |>
        dplyr::mutate(
            mean_degree = mean(centrality_degree, na.rm = TRUE),
            mean_eigen = mean(centrality_eigen, na.rm = TRUE),
            mean_closeness = mean(centrality_closeness, na.rm = TRUE),
            mean_betweenness = mean(centrality_betweenness, na.rm = TRUE),
            n_features = dplyr::n()
        ) |>
        dplyr::ungroup()

    result <- list(
        exposure_impact = impact_df
    )

    # Store or return
    if (action == "add") {
        all_metadata <- MultiAssayExperiment::metadata(expomicset)
        all_metadata$network$exposure_impact[[feature_type]] <- result
        MultiAssayExperiment::metadata(expomicset) <- all_metadata

        step_record <- list(
            run_exposure_impact = list(
                timestamp = Sys.time(),
                params = list(feature_type = feature_type),
                notes = paste(
                    "Computed exposure impact using",
                    feature_type,
                    "correlation network."
                )
            )
        )

        MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(expomicset)$summary$steps,
            step_record
        )
        return(expomicset)
    } else {
        return(result)
    }
}
