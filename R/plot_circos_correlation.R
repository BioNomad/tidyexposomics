#' Plot Circular Network of Exposure Relationships
#'
#' Generates a circular network plot to visualize relationships between
#' exposures, either based on correlation ("exposures") or
#' shared features ("degs", "factors").
#'
#' @param exposomicset A MultiAssayExperiment object.
#' @param feature_type One of "exposures", "degs", or "factors".
#' @param exposure_cols Character vector of exposures to include
#' (only for "exposures").
#' @param corr_threshold Minimum |correlation| (only for "exposures").
#' @param shared_cutoff Minimum number of shared features
#' (only for "degs" or "factors"). Default = 10.
#' @param annotation_colors Optional named vector of colors for categories.
#' @param low low value color for edges.
#' @param mid middle value color for edges.
#' @param high high value color for edges.
#' @param midpoint Midpoint for edge color gradient. Defaults to 0
#' (for correlations) or mean shared features.
#'
#' @return A ggplot object (ggraph circular plot).
#'
#' @import ggplot2
#' @importFrom dplyr filter rename arrange
#' @importFrom igraph graph_from_data_frame
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run correlation analysis
#' mae <- mae |>
#'     run_correlation(
#'         feature_type = "exposures",
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' # create the circos plot
#' circos_plot <- mae |>
#'     plot_circos_correlation(
#'         feature_type = "exposures"
#'     )
#'
#' @export
plot_circos_correlation <- function(
  exposomicset,
  feature_type = c(
      "degs",
      "omics",
      "factors",
      "factor_features",
      "exposures", "pcs"
  ),
  exposure_cols = NULL,
  corr_threshold = NULL,
  shared_cutoff = 10,
  annotation_colors = NULL,
  low = "#006666",
  mid = "white",
  high = "#8E0152",
  midpoint = NULL
) {
    # require(ggplot2)
    # require(ggraph)
    # require(igraph)
    .check_suggested(pkg = "ggraph")
    .check_suggested(pkg = "scales")

    feature_type <- match.arg(feature_type)

    # Load correlation or shared feature data
    if (feature_type == "exposures") {
        correlation_df <- MultiAssayExperiment::metadata(exposomicset) |>
            purrr::pluck("correlation", "exposures")

        if (is.null(correlation_df)) {
            stop("Run `run_correlation(..., feature_type = 'exposures')` first.")
        }

        # Optional exposure filtering
        if (!is.null(exposure_cols)) {
            correlation_df <- dplyr::filter(
                correlation_df,
                var1 %in% exposure_cols & var2 %in% exposure_cols
            )
        }

        if (!is.null(corr_threshold)) {
            correlation_df <- dplyr::filter(
                correlation_df, abs(correlation) >= corr_threshold
            )
        }

        correlation_df$edge_weight <- abs(correlation_df$correlation)
        correlation_df <- dplyr::rename(correlation_df,
            source = var1,
            target = var2
        )
        edge_color_var <- correlation_df$correlation
        midpoint_val <- ifelse(is.null(midpoint), 0, midpoint)
    } else {
        # degs or factors
        tag <- if (feature_type == "degs") {
            "omics_exposure_deg_correlation"
        } else {
            "omics_exposure_factor_correlation"
        }

        correlation_list <- MultiAssayExperiment::metadata(exposomicset)$correlation
        correlation_df <- correlation_list[[feature_type]]

        # if (is.null(correlation_df)){
        #   stop(paste0(
        #     "No correlation found for feature_type = '",
        #     feature_type,
        #     "'. Run `run_correlation(..., feature_type = '",
        #     feature_type,
        #     "')` first."))
        # }
        if (is.null(correlation_df)) {
            stop(sprintf(
                "No correlation found for feature_type = '%s'.
        Run `run_correlation(..., feature_type = '%s')` first.",
                feature_type, feature_type
            ))
        }

        feature_list <- split(
            correlation_df,
            correlation_df$exposure
        ) |>
            purrr::map(~ unique(.x$feature))

        overlap_df <- .get_pairwise_overlaps(feature_list)
        overlap_df <- overlap_df |>
            dplyr::filter(num_shared > shared_cutoff)

        correlation_df <- overlap_df |>
            dplyr::rename(
                source = source,
                target = target,
                edge_weight = num_shared
            )

        edge_color_var <- correlation_df$edge_weight
        midpoint_val <- ifelse(is.null(midpoint),
            mean(correlation_df$edge_weight),
            midpoint
        )
    }

    # Get node metadata
    codebook <- MultiAssayExperiment::metadata(exposomicset)$codebook
    node_vars <- unique(c(correlation_df$source, correlation_df$target))
    vertex_df <- dplyr::filter(codebook, variable %in% node_vars)

    # Sort Correlation Df
    correlation_df <- correlation_df |>
        arrange(abs(edge_weight))

    # Ensure igraph object
    graph <- igraph::graph_from_data_frame(correlation_df,
        vertices = vertex_df,
        directed = FALSE
    )

    # Colors
    all_categories <- na.omit(unique(vertex_df$category))
    if (!is.null(annotation_colors)) {
        if (!all(all_categories %in% names(annotation_colors))) {
            stop("All node categories must be included in annotation_colors.")
        }
        cat_colors <- annotation_colors
    } else {
        cat_colors <- setNames(
            tidy_exp_pal[seq_along(all_categories)],
            all_categories
        )
    }

    # Plot
    ggraph::ggraph(graph,
        layout = "linear",
        circular = TRUE
    ) +
        ggraph::geom_edge_arc(aes(
            color = edge_weight,
            width = scales::rescale(edge_weight, to = c(1, 5))
        )) +
        ggraph::geom_node_point(aes(fill = category),
            shape = 21,
            size = 4,
            alpha = 0.8
        ) +
        ggraph::geom_node_text(
            aes(
                label = name,
                x = 1.1 * x,
                y = 1.1 * y,
                angle = ggraph::node_angle(x, y)
            ),
            hjust = "outward",
            fontface = "bold.italic",
            check_overlap = TRUE
        ) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::scale_fill_manual(values = cat_colors) +
        ggraph::scale_edge_color_gradient2(
            low = low,
            mid = mid,
            high = high,
            midpoint = midpoint_val,
            guide = ggraph::guide_edge_colorbar()
        ) +
        ggplot2::coord_fixed(xlim = c(-2, 2), ylim = c(-2, 2)) +
        ggplot2::guides(edge_width = "none") +
        ggplot2::labs(
            edge_color = ifelse(feature_type == "exposures",
                "Correlation",
                "Shared Features"
            ),
            fill = "Category"
        )
}
