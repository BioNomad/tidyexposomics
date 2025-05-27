#' Plot Overlap of Top Factor Features Across Experiments
#'
#' Visualizes the **overlap of top features** across multi-omics integration factors using barplots and a Venn diagram.
#'
#' @param expomicset A `MultiAssayExperiment` object with `top_factor_features` and `common_top_factor_features` in metadata.
#' @param venn_text_size Text size for Venn diagram labels. Default is `3`.
#' @param venn_stroke_size Outline stroke width in Venn diagram. Default is `0.1`.
#' @param venn_set_name_size Text size for set names in the Venn diagram. Default is `3.5`.
#' @param venn_show_percent Logical, whether to show percentages in the Venn diagram. Default is `FALSE`.
#' @param venn_colors Optional named vector of colors for each factor in the Venn diagram.
#' @param shared_bar_colors Optional named vector of colors for the shared/unique barplot fill.
#' @param da_bar_colors Optional named vector of colors for DEG (differentially abundant) barplots.
#' @param da_bar_facet_cols Optional vector of colors for facet strip backgrounds in DEG-related barplots.
#' @param da_bar_facet_alpha Numeric from 0–1 for facet background transparency. Default is `0.5`.
#' @param left_heights A numeric vector of relative heights for the Venn + shared DEG barplot panel. Default is `c(2, 1)`.
#' @param col_widths A numeric vector of relative widths for left (Venn/shared) vs right (DEG barplot) columns. Default is `c(1, 1)`.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Displays overlap of top features across integration factors using a Venn diagram (`ggvenn`).
#'   \item Visualizes barplots of:
#'     \itemize{
#'       \item Shared vs unique top features across experiments.
#'       \item Differentially abundant top features (per factor).
#'       \item DEGs among only the shared features (overall).
#'     }
#'   \item Adds optional customization for layout spacing, colors, and transparency.
#' }
#'
#' @return A `patchwork` object showing the combined plots: barplots and Venn diagram.
#'
#' @examples
#' \dontrun{
#' plot_factor_overlap(expomicset, venn_show_percent = TRUE)
#' }
#'
#' @export

plot_factor_overlap <- function(
    expomicset,
    venn_text_size = 3,
    venn_stroke_size=0.1,
    venn_set_name_size = 3.5,
    venn_show_percent=FALSE,
    venn_colors = NULL,
    shared_bar_colors = NULL,
    da_bar_colors = NULL,
    da_bar_facet_cols = NULL,
    da_bar_facet_alpha=0.5,
    left_heights = c(2, 1),
    col_widths = c(1, 1)
){
  # Check if top factor features are available
  if(!("top_factor_features" %in% names(MultiAssayExperiment::metadata(expomicset)))){
    stop("Please run `extract_top_factor_features()` first to generate the top factor features.")
  }

  # Extract top factor features
  top_factor_features <- MultiAssayExperiment::metadata(expomicset)$top_factor_features

  # Extract common feature information
  common_features <- MultiAssayExperiment::metadata(expomicset)$common_top_factor_features

  if(!is.null(venn_colors)){
    venn_colors <- venn_colors
  } else{
    venn_colors <- get_palette("npg", k = length(unique(top_factor_features$factor)))
  }
  # Venn Diagram of top factor features
  venn <- top_factor_features |>
    (\(df)split(df,df$factor))() |>
    map(~.x |>
          pull(feature) |>
          unique()) |>
    ggvenn::ggvenn(text_size = venn_text_size,
                   stroke_size = venn_stroke_size,
                   set_name_size = venn_set_name_size,
                   show_percent = venn_show_percent)+
    scale_fill_manual(values=venn_colors)+
    theme_void() +  # remove background, grid, etc.
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank(),
      aspect.ratio = 1

    )+
    coord_cartesian(clip = "off")

  if(!is.null(shared_bar_colors)){
    shared_bar_colors <- shared_bar_colors
  } else{
    shared_bar_colors <- get_palette("lancet", k = length(unique(top_factor_features$exp_name)))
  }


  if(!is.null(da_bar_colors)){
    da_bar_colors <- da_bar_colors
  } else{
    da_bar_colors <- get_palette("uchicago", k = length(unique(top_factor_features$exp_name)))
  }

  if(!is.null(da_bar_facet_cols)){
    da_bar_facet_cols <- da_bar_facet_cols
  } else{
    # Use a default palette for facets if not provided
    da_bar_facet_cols <- get_palette("uchicago", k = length(unique(top_factor_features$factor)))
  }

  # Barplot of experiment names, and colored by if the feature is shared
  barplot_shared <- top_factor_features |>
    dplyr::select(exp_name,factor,feature)|>
    distinct() |>
    mutate(exp_name_feature = paste0(exp_name, "_", feature)) |>
    dplyr::mutate(is_shared = ifelse(
      exp_name_feature %in%
        common_features$exp_name_feature,
      "Shared",
      "Unique")) |>
    group_by(exp_name, factor,is_shared) |>
    summarise(n = n()) |>
    mutate(exp_name = paste0(
      exp_name, " (",
      sum(n[is_shared == "Shared"], na.rm = TRUE), "/",
      sum(n),
      ")"
    )) |>
    ungroup() |>
    ggplot(aes(
      x=n,
      y=forcats::fct_reorder(exp_name, n),
      fill=is_shared
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    ggh4x::facet_grid2(
      factor~.,
      scales = "free_y",
      space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = scales::alpha(
            da_bar_facet_cols,
            da_bar_facet_alpha
          )
        )
      )
    ) +
    scale_fill_manual(values=shared_bar_colors)+
    ggpubr::theme_pubr(legend="right",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"))+
    guides(color="none")+
    labs(title = "Top Factor Feature Overlap",
         y = "",
         x = "No. of Features",
         fill = ""
    )

  # Barplot of experiment names, and colored by if the feature is deg
  barplot_deg <- top_factor_features |>
    dplyr::select(exp_name,factor,feature)|>
    distinct() |>
    mutate(exp_name_feature = paste0(exp_name, "_", feature)) |>
    dplyr::mutate(is_deg = ifelse(
      exp_name_feature %in%
        common_features$exp_name_feature[
          common_features$is_deg == TRUE],
      "Differentially Abundant",
      "Not Differentially Abundant")) |>
    group_by(exp_name, factor,is_deg) |>
    summarise(n = n()) |>
    mutate(exp_name = paste0(
      exp_name, " (",
      sum(n[is_deg == "Differentially Abundant"], na.rm = TRUE), "/",
      sum(n),
      ")"
    )) |>
    ungroup() |>
    ggplot(aes(
      x=n,
      y=forcats::fct_reorder(exp_name, n),
      fill=is_deg
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    ggh4x::facet_grid2(
      factor~.,
      scales = "free_y",
      space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = scales::alpha(
            da_bar_facet_cols,
            da_bar_facet_alpha
          )
        )
      )
    ) +
    scale_fill_manual(values=da_bar_colors)+
    ggpubr::theme_pubr(legend="right",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"),
          strip.text.y = element_text(angle=0,face = "bold.italic"))+
    guides(color="none")+
    labs(title = "Differentially Abundant Top Factor Features",
         y = "",
         x = "No. of Features",
         fill = ""
    )

  # Barplot of shared features that are also DEGs (bottom-right)
  barplot_shared_deg <- common_features |>
    dplyr::select(exp_name,feature)|>
    distinct() |>
    mutate(exp_name_feature = paste0(exp_name, "_", feature)) |>
    dplyr::mutate(is_deg = ifelse(
      exp_name_feature %in%
        common_features$exp_name_feature[
          common_features$is_deg == TRUE],
      "Differentially Abundant",
      "Not Differentially Abundant")) |>
    group_by(exp_name,is_deg) |>
    summarise(n = n()) |>
    mutate(exp_name = paste0(
      exp_name, " (",
      sum(n[is_deg == "Differentially Abundant"], na.rm = TRUE), "/",
      sum(n),
      ")"
    )) |>
    ungroup() |>
    ggplot(aes(
      x=n,
      y=forcats::fct_reorder(exp_name, n),
      fill=is_deg
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    scale_fill_manual(values=da_bar_colors)+
    ggpubr::theme_pubr(legend="right",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"),
          strip.text.y = element_text(angle=0,face = "bold.italic"))+
    guides(color="none")+
    labs(title = "Differentially Abundant Shared Features",
         y = "",
         x = "No. of Features",
         fill = ""
    )


  # Top row: barplot_shared and barplot_deg (equal heights)
  right_p <- barplot_shared | barplot_deg
  right_p <-  barplot_deg


  # Bottom row: Venn and the shared DEG barplot (6:1 height ratio)
  left_p <- wrap_elements(full = venn) / barplot_shared_deg +
    plot_layout(heights = left_heights)

  # Combine into full plot
  final_plot <- (left_p | right_p)+
    plot_layout(widths = col_widths)

  final_plot
}

