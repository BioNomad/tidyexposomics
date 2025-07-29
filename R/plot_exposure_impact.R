#' Plot Exposure Impact on Network Centrality Metrics
#'
#' Visualizes the impact of exposures on network centrality measures of associated features
#' (e.g., genes or latent factors) as a heatmap. Each exposure is scored by four centrality metrics,
#' scaled within metric, and grouped by exposure category.
#'
#' @param expomicset A `MultiAssayExperiment` object with results from `run_exposure_impact()`.
#' @param feature_type Character string specifying the feature type. One of `"degs"`, `"omics"`, or `"factors"`.
#' @param min_per_group Minimum number of features per exposure for inclusion (not currently used). Default is `5`.
#' @param facet_cols Optional named vector of colors for exposure categories.
#' @param bar_cols Optional vector of colors for bar plots (if enabled).
#' @param alpha Transparency level for category strips (if enabled). Default is `0.3`.
#' @param ncol,nrow Layout for optional patchwork combination (currently unused). Default: `ncol = 2`, `nrow = 1`.
#' @param heights, widths Relative heights and widths for combined plots (currently unused). Defaults: `c(1,1)`, `c(2,1)`.
#'
#' @return A `ggplot`/`patchwork` object showing a heatmap of scaled network centrality scores per exposure, annotated by category.
#'
#' @details
#' This function uses the output of `run_exposure_impact()` to calculate and visualize the mean centrality
#' values for each exposure across its associated features. The following network centrality metrics are shown:
#' \itemize{
#'   \item Degree centrality
#'   \item Eigenvector centrality
#'   \item Closeness centrality
#'   \item Betweenness centrality
#' }
#' All values are scaled within metric across exposures. A side bar indicates the category of each exposure.
#'
#' @examples
#' \dontrun{
#' plot_exposure_impact(expomicset, feature_type = "degs")
#' }
#'
#' @export
plot_exposure_impact <- function(
    expomicset,
    feature_type = c("degs", "omics", "factors"),
    min_per_group=5,
    facet_cols=NULL,
    bar_cols=NULL,
    alpha=0.3,
    ncol=2,
    nrow=1,
    heights = c(1,1),
    widths = c(2,1)){

  require(ggplot2)
  require(patchwork)

  # Check that exposure impact analysis has been run
  if(!("exposure_impact" %in% names(MultiAssayExperiment::metadata(expomicset)$network))){
    stop("Please run `run_exposure_impact()` first.")
  }

  exposure_impact_degree <-
    expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("network") |>
    purrr::pluck("exposure_impact") |>
    purrr::pluck(feature_type) |>
    purrr::pluck("exposure_impact")

  # exposure_impact_deg_n <-
  #   expomicset |>
  #   MultiAssayExperiment::metadata() |>
  #   purrr::pluck("network") |>
  #   purrr::pluck("exposure_impact") |>
  #   purrr::pluck(feature_type) |>
  #   purrr::pluck("deg_association_summary")

  # set facet colors
  if(!is.null(facet_cols)){
    facet_cols <- facet_cols
  } else{
    facet_cols <- tidy_exp_pal[1:length(
      unique(
        exposure_impact_degree |>
          dplyr::pull(category)
      ))]
  }

  # set bar colors
  if(!is.null(bar_cols)){
    bar_cols <- bar_cols
  } else{
    bar_cols <- ggsci::pal_uchicago("default")(length(
      unique(
        exposure_impact_degree |>
          dplyr::pull(exp_name)
      )))
  }

  # Create a boxplot of the degree centrality of genes that are associated
  # with a particular exposure
  # boxplot <- exposure_impact_degree |>
  #   dplyr::filter(n_features >min_per_group) |>
  #   ggplot(aes(
  #     x=reorder(exposure,-mean_degree),
  #     y=centrality_degree,
  #     fill=category,
  #     color=category
  #   ))+
  #   geom_boxplot(alpha=0.7) +
  #   geom_jitter(alpha=0.1)+
  #   ggpubr::rotate_x_text(angle=45)+
  #   theme_minimal()+
  #   #ggpubr::theme_pubclean()+
  #   scale_fill_manual(values = facet_cols) +
  #   ggpubr::rotate_x_text(angle=45)+
  #   ggh4x::facet_grid2(
  #     ~category,
  #     scales = "free_x",
  #     space = "free_x",
  #     strip = ggh4x::strip_themed(
  #       background_x = ggh4x::elem_list_rect(
  #         fill = scales::alpha(
  #           facet_cols,
  #           alpha
  #         )
  #       )
  #     )
  #   ) +
  #   scale_color_manual(values=facet_cols)+
  #   theme(strip.text.x = element_text(face = "bold.italic"),
  #         legend.position = "none")+
  #   labs(
  #     x="",
  #     y="Degree",
  #   )

  # barplot <- exposure_impact_deg_n |>
  #   dplyr::select(-percent_assoc_exp) |>
  #   dplyr::mutate(left_over= total_deg-n_assoc_exp) |>
  #   dplyr::select(-total_deg) |>
  #   tidyr::pivot_longer(-exp_name,
  #                       names_to = "group",
  #                       values_to = "value") |>
  #   dplyr::mutate(group=dplyr::case_when(
  #     group == "n_assoc_exp" ~ "Associated with Exposures",
  #     group == "left_over" ~ "Not Associated with Exposures"
  #   )) |>
  #   ggplot(aes(
  #     x=value,
  #     y=reorder(exp_name,value),
  #     fill=group
  #   ))+
  #   geom_bar(stat="identity")+
  #   ggpubr::rotate_x_text(angle=45)+
  #   ggpubr::theme_pubr(legend="right")+
  #   scale_fill_manual(values = bar_cols)+
  #   labs(
  #     x="Number of Features",
  #     y="",
  #     fill=""
  #   )

  centrality_heatmap_df <- exposure_impact_degree |>
    dplyr::select(exposure,
                  category,
                  mean_degree,
                  mean_eigen,
                  mean_closeness,
                  mean_betweenness) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(
      cols = -c(exposure,category),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::group_by(metric) |>
    dplyr::mutate(scaled_value = scale(value)[, 1]) |>  # scale() returns matrix
    dplyr::ungroup()

  # Step 1: create exposure order based on centrality_heatmap_df
  exposure_order <- centrality_heatmap_df |>
    dplyr::group_by(exposure) |>
    dplyr::summarise(avg = mean(scaled_value, na.rm = TRUE)) |>
    dplyr::arrange(avg) |>
    dplyr::pull(exposure)

  # Step 2: apply this order in both plots
  centrality_heatmap_df <- centrality_heatmap_df |>
    dplyr::mutate(exposure = factor(exposure, levels = exposure_order))

  category_bar <- centrality_heatmap_df |>
    dplyr::select(exposure, category) |>
    dplyr::distinct() |>
    ggplot(aes(y = exposure, x = 1, fill = category)) +
    geom_tile() +
    scale_fill_manual(values = facet_cols, name = "Category") +
    theme_void() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10)
    )

  heatmap_plot <- centrality_heatmap_df |>
    mutate(metric = dplyr::case_when(
      metric == "mean_degree" ~ "Mean Degree",
      metric == "mean_eigen" ~ "Mean Eigenvector Centrality",
      metric == "mean_closeness" ~ "Mean Closeness Centrality",
      metric == "mean_betweenness" ~ "Mean Betweenness Centrality"
    )) |>
    ggplot(aes(x = metric,
               y = exposure,
               fill = scaled_value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue4",
      mid = "white",
      high = "red4",
      midpoint = 0,
      name = "Scaled Value"
    ) +
    labs(x = "", y = "", fill = "Scaled\nValue") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(2, 0, 2, 2)
    )
  heatmap_plot <- category_bar + heatmap_plot +
    plot_layout(ncol = 2, widths = c(0.15, 1),guides = "collect")

  # heatmap_plot <- centrality_heatmap_df |>
  #   mutate(metric = case_when(
  #     metric == "mean_degree" ~ "Mean Degree",
  #     metric == "mean_eigen" ~ "Mean Eigenvector Centrality",
  #     metric == "mean_closeness" ~ "Mean Closeness Centrality",
  #     metric == "mean_betweenness" ~ "Mean Betweenness Centrality"
  #   )) |>
  #   ggplot(aes(x = metric,
  #              y = reorder(exposure, -scaled_value),
  #              fill = scaled_value)) +
  #   geom_tile(color = "white") +
  #   ggh4x::facet_grid2(
  #     category~.,
  #     scales = "free_y",
  #     space = "free_y",
  #     strip = ggh4x::strip_themed(
  #       background_x = ggh4x::elem_list_rect(
  #         fill = scales::alpha(
  #           facet_cols,
  #           alpha
  #         )
  #       )
  #     )
  #   ) +
  #   scale_fill_gradient2(
  #     low = "blue",
  #     mid = "white",
  #     high = "red",
  #     midpoint = 0,
  #     name = "Scaled\nValue"
  #   ) +
  #   labs(x = "",
  #        y = "",
  #        fill = "Scaled\nValue") +
  #   theme_minimal(base_size = 11) +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     axis.text.y = element_text(face = "bold"),
  #     panel.grid = element_blank(),
  #     strip.text.y = element_text(angle=0)
  #   )

  # final_plot <- patchwork::wrap_plots(
  #   boxplot,
  #   heatmap_plot,
  #   ncol = ncol,
  #   nrow = nrow)+
  #   patchwork::plot_layout(
  #     heights = heights,
  #     widths = widths
  #   )


  return(heatmap_plot)
}
