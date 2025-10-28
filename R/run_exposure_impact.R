#' Calculate Exposure Impact from Feature-Exposure Correlation Networks
#'
#' Generalized centrality-based exposure impact analysis using DEG, omics,
#' or factor features.
#'
#' @param exposomicset A `MultiAssayExperiment` object with correlation
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
  exposomicset,
  feature_type = c(
      "degs",
      "omics",
      "factor_features"
  ),
  pval_col = "adj.P.Val",
  pval_thresh = 0.1,
  action = c("add", "get")
) {
    .check_suggested(pkg = "tidygraph")

    feature_type <- match.arg(feature_type)
    action <- match.arg(action)

    md <- MultiAssayExperiment::metadata(exposomicset)

    # --- checks
    if (is.null(md$correlation[[feature_type]])) {
        stop(
            "Correlation results missing for feature_type = ", feature_type,
            ". Please run run_correlation(..., feature_type = '", feature_type, "')."
        )
    }

    # prefer "..._feature_cor" but fall back to plain "network_<feature_type>"
    graph_key <- paste0("network_", feature_type, "_feature_cor")
    if (is.null(md$network[[graph_key]])) {
        fallback_key <- paste0("network_", feature_type)
        if (!is.null(md$network[[fallback_key]])) {
            graph_key <- fallback_key
        } else {
            stop(
                "Network missing for feature_type = ", feature_type,
                ". Please run run_create_network(..., feature_type = '", feature_type, "')."
            )
        }
    }

    correlation_df <- md$correlation[[feature_type]]
    net_graph <- md$network[[graph_key]]$graph

    # Optional p-val filter if column exists
    if (!is.null(pval_col) && pval_col %in% names(correlation_df)) {
        correlation_df <- dplyr::filter(correlation_df, .data[[pval_col]] <= pval_thresh)
    }

    # --- centralities from graph
    node_centrality <- tidygraph::as_tbl_graph(net_graph) |>
        tidygraph::activate(nodes) |>
        dplyr::mutate(
            centrality_degree      = tidygraph::centrality_degree(),
            centrality_eigen       = tidygraph::centrality_eigen(),
            centrality_closeness   = tidygraph::centrality_closeness(),
            centrality_betweenness = tidygraph::centrality_betweenness()
        ) |>
        tibble::as_tibble()

    # --- construct join keys that mirror create_network()
    # Case A: exposure-feature schema: feature, exposure, exp_name, category
    if (all(c("feature", "exposure", "exp_name", "category") %in% names(correlation_df))) {
        df <- correlation_df |>
            dplyr::mutate(
                feature_vertex_id  = paste(.data$feature, .data$exp_name, sep = "_"),
                exposure_vertex_id = paste(.data$exposure, .data$category, sep = "_")
            )

        exposure_col <- "exposure"

        # generic var1/var2 with explicit types (expect var1 == exposure)
    } else if (all(c("var1", "var2", "var1_type", "var2_type") %in% names(correlation_df))) {
        df <- correlation_df |>
            dplyr::filter(.data$var1_type == "exposure", .data$var2_type != "exposure") |>
            dplyr::mutate(
                exposure            = .data$var1,
                feature             = .data$var2,
                feature_vertex_id   = paste(.data$var2, .data$var2_type, sep = "_"),
                exposure_vertex_id  = paste(.data$var1, .data$var1_type, sep = "_")
            )

        exposure_col <- "exposure"

        # Case C: cross-experiment pairings (use exp_name_1/exp_name_2 as groups)
    } else if (all(c("var1", "var2", "exp_name_1", "exp_name_2") %in% names(correlation_df))) {
        df <- correlation_df |>
            dplyr::mutate(
                # Naming here keeps "exposure" just as a grouping label for summary;
                # adapt as needed for your schema.
                exposure            = .data$var1,
                feature             = .data$var2,
                feature_vertex_id   = paste(.data$var2, .data$exp_name_2, sep = "_"),
                exposure_vertex_id  = paste(.data$var1, .data$exp_name_1, sep = "_")
            )

        exposure_col <- "exposure"
    } else {
        stop("Unrecognized correlation format: cannot construct vertex ids to match the graph.")
    }

    # --- join feature nodes to their centralities via vertex_id
    df2 <- df |>
        dplyr::left_join(
            node_centrality,
            by = c("feature_vertex_id" = "name") # name is vertex_id in the graph
        )

    # --- aggregate to exposure-level impact
    impact_df <- df2 |>
        dplyr::filter(!is.na(.data$centrality_degree)) |>
        dplyr::group_by(.data[[exposure_col]]) |>
        dplyr::summarise(
            mean_degree = mean(centrality_degree, na.rm = TRUE),
            mean_eigen = mean(centrality_eigen, na.rm = TRUE),
            mean_closeness = mean(centrality_closeness, na.rm = TRUE),
            mean_betweenness = mean(centrality_betweenness, na.rm = TRUE),
            n_features = dplyr::n(),
            .groups = "drop"
        ) |>
        left_join(
            exposomicset@metadata$codebook,
            by = c("exposure" = "variable")
        )

    result <- list(exposure_impact = impact_df)

    if (action == "add") {
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$network$exposure_impact[[feature_type]] <- result
        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

        step_record <- list(run_exposure_impact = list(
            timestamp = Sys.time(),
            params    = list(feature_type = feature_type, pval_col = pval_col, pval_thresh = pval_thresh),
            notes     = paste("Computed exposure impact using", feature_type, "correlation network.")
        ))
        MultiAssayExperiment::metadata(exposomicset)$summary$steps <-
            c(MultiAssayExperiment::metadata(exposomicset)$summary$steps, step_record)

        return(exposomicset)
    } else {
        return(result)
    }
}
