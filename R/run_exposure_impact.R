#' Calculate Exposure Impact from Feature-Exposure Correlation Networks
#'
#' Generalized centrality-based exposure impact analysis using DEG, omics, or factor features.
#'
#' @param expomicset A `MultiAssayExperiment` object with correlation and network metadata.
#' @param feature_type One of `"degs"`, `"omics"`, or `"factor_features"`.
#' @param pval_col Column in differential abundance results to filter DEGs. Default = `"adj.P.Val"`.
#' @param pval_thresh DEG p-value threshold. Ignored unless `feature_type == "degs"`.
#' @param action Either `"add"` (store in metadata) or `"get"` (return list).
#'
#' @return Either an updated MultiAssayExperiment (if action = "add") or a list.
#' @export
run_exposure_impact <- function(
    expomicset,
    feature_type = c("degs",
                     "omics",
                     "factor_features"),
    pval_col = "adj.P.Val",
    pval_thresh = 0.1,
    action = c("add", "get")
) {
  feature_type <- match.arg(feature_type)
  action <- match.arg(action)
  md <- MultiAssayExperiment::metadata(expomicset)

  # Validate required correlation and network
  if (is.null(md$correlation[[feature_type]])) {
    stop("Correlation results missing for feature_type = ", feature_type,
         ". Please run run_correlation(..., feature_type = '", feature_type, "').")
  }
  if (is.null(md$network[[paste0("network_", feature_type)]])) {
    stop("Network missing for feature_type = ", feature_type,
         ". Please run run_create_network(..., feature_type = '", feature_type, "').")
  }

  # Get correlation and network
  correlation_df <- md$correlation[[feature_type]]
  net_graph <- md$network[[paste0("network_", feature_type,"_feature_cor")]]$graph

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
      dplyr::rename(feature = feature,
                    exposure = exposure)

  } else if ("var1" %in% names(correlation_df)) {
    # Handle generalized format
    correlation_df <- correlation_df |>
      dplyr::filter(var1_type == "exposure",
                    var2_type != "exposure") |>
      dplyr::rename(exposure = var1,
                    feature = var2)
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
    MultiAssayExperiment::metadata(expomicset)$network$exposure_impact[[feature_type]] <- result

    step_record <- list(
      run_exposure_impact = list(
        timestamp = Sys.time(),
        params = list(feature_type = feature_type),
        notes = paste("Computed exposure impact using",
                      feature_type,
                      "correlation network.")
      ))

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )
    return(expomicset)
  } else {
    return(result)
  }
}
# run_exposure_impact <- function(
#     expomicset,
#     feature_type = c("degs",
#                      "omics",
#                      "factors",
#                      "factor_features",
#                      "exposures",
#                      "degs_feature_cor",
#                      "omics_feature_cor",
#                      "factor_features_feature_cor"),
#     pval_col = "adj.P.Val",
#     pval_thresh = 0.1,
#     action = c("add", "get")
# ) {
#   feature_type <- match.arg(feature_type)
#   action <- match.arg(action)
#   md <- MultiAssayExperiment::metadata(expomicset)
#
#   # Validate required correlation and network
#   if (is.null(md$correlation[[feature_type]])) {
#     stop("Correlation results missing for feature_type = ", feature_type,
#          ". Please run run_correlation(..., feature_type = '", feature_type, "').")
#   }
#   if (is.null(md$network[[paste0("network_", feature_type)]])) {
#     stop("Network missing for feature_type = ", feature_type,
#          ". Please run run_create_network(..., feature_type = '", feature_type, "').")
#   }
#
#   # Get correlation and network
#   correlation_df <- md$correlation[[feature_type]]
#   net_graph <- md$network[[paste0("network_", feature_type)]]$graph
#
#   # Compute centrality for all nodes
#   # node_centrality <- net_graph |>
#   #   tidygraph::as_tbl_graph() |>
#   #   tidygraph::activate(nodes) |>
#   #   dplyr::mutate(centrality = tidygraph::centrality_degree()) |>
#   #   dplyr::as_tibble()
#
#   node_centrality <- net_graph |>
#     tidygraph::as_tbl_graph() |>
#     tidygraph::activate(nodes) |>
#     dplyr::mutate(
#       centrality_degree = tidygraph::centrality_degree(),
#       centrality_eigen = tidygraph::centrality_eigen(),
#       centrality_closeness = tidygraph::centrality_closeness(),
#       centrality_betweenness = tidygraph::centrality_betweenness()
#     ) |>
#     dplyr::as_tibble()
#
#   # Merge with correlation edges (exposure-feature)
#   if (feature_type == "degs") {
#     correlation_df <- correlation_df |>
#       dplyr::rename(feature = feature,
#                     exposure = exposure)
#
#   } else if ("var1" %in% names(correlation_df)) {
#     # Handle generalized format
#     correlation_df <- correlation_df |>
#       dplyr::filter(var1_type == "exposure",
#                     var2_type != "exposure") |>
#       dplyr::rename(exposure = var1,
#                     feature = var2)
#   }
#
#   # impact_df <- correlation_df |>
#   #   dplyr::left_join(node_centrality,
#   #                    by = c("feature" = "name")) |>
#   #   dplyr::filter(!is.na(centrality)) |>
#   #   dplyr::group_by(exposure) |>
#   #   dplyr::mutate(
#   #     mean_centrality = mean(centrality, na.rm = TRUE),
#   #     n_features = dplyr::n()
#   #   ) |>
#   #   dplyr::ungroup()
#
#   impact_df <- correlation_df |>
#     dplyr::left_join(node_centrality, by = c("feature" = "name")) |>
#     dplyr::filter(!is.na(centrality_degree)) |>
#     dplyr::group_by(exposure) |>
#     dplyr::mutate(
#       mean_degree = mean(centrality_degree, na.rm = TRUE),
#       mean_eigen = mean(centrality_eigen, na.rm = TRUE),
#       mean_closeness = mean(centrality_closeness, na.rm = TRUE),
#       mean_betweenness = mean(centrality_betweenness, na.rm = TRUE),
#       n_features = dplyr::n()
#     ) |>
#     dplyr::ungroup()
#
#   # Optionally compute DEG counts per exposure
#   if (feature_type == "degs") {
#     if (is.null(md$differential_analysis$differential_abundance)) {
#       stop("Differential abundance results missing for DEG impact summary.")
#     }
#
#     da_df <- md$differential_analysis$differential_abundance
#
#     deg_summary <- da_df |>
#       dplyr::filter(!!rlang::sym(pval_col) < pval_thresh) |>
#       dplyr::left_join(
#         correlation_df |>
#           dplyr::select(feature, exp_name) |>
#           dplyr::distinct() |>
#           dplyr::mutate(exposure_association = "yes"),
#         by = c("feature", "exp_name")
#       ) |>
#       dplyr::mutate(
#         exposure_association = ifelse(is.na(exposure_association), "no", exposure_association)
#       ) |>
#       dplyr::group_by(exp_name) |>
#       dplyr::summarise(
#         total_deg = dplyr::n(),
#         n_assoc_exp = sum(exposure_association == "yes"),
#         percent_assoc_exp = n_assoc_exp / total_deg * 100
#       )
#
#     result <- list(
#       exposure_impact = impact_df,
#       deg_association_summary = deg_summary
#     )
#   } else {
#     result <- list(
#       exposure_impact = impact_df
#     )
#   }
#
#   # Store or return
#   if (action == "add") {
#     MultiAssayExperiment::metadata(expomicset)$network$exposure_impact[[feature_type]] <- result
#
#     step_record <- list(
#       run_exposure_impact = list(
#       timestamp = Sys.time(),
#       params = list(feature_type = feature_type),
#       notes = paste("Computed exposure impact using",
#                     feature_type,
#                     "correlation network.")
#     ))
#
#     MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
#       MultiAssayExperiment::metadata(expomicset)$summary$steps,
#       step_record
#     )
#     return(expomicset)
#   } else {
#     return(result)
#   }
# }
