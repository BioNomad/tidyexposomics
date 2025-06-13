#' Correlation Heatmap of Exposure Variables with Category Annotations
#'
#' Generates an upper-triangle **heatmap** of correlations between exposure variables,
#' with side annotation bars indicating variable categories.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure data and correlation results.
#' @param exposure_cols A character vector of exposure variable names to include in the plot. If `NULL`, includes all.
#' @param corr_threshold A numeric value for filtering correlations by absolute value (e.g., `0.2`).
#' @param annotation_colors Optional named character vector of colors for exposure categories.
#' @param low Color for strong negative correlations. Default is `"#006666"`.
#' @param mid Color for zero correlation. Default is `"white"`.
#' @param high Color for strong positive correlations. Default is `"#8E0152"`.
#' @param midpoint Midpoint for the color scale. Default is `0`.
#'
#' @details
#' This function:
#' - Uses precomputed correlations from `metadata(expomicset)$exposure_correlation$correlation_table`.
#' - Filters to an upper triangle matrix (i.e., only combinations where `var1 <= var2`).
#' - Applies optional thresholding with `corr_threshold` to simplify visualization.
#' - Adds **top and left annotation bars** colored by `category_1` and `category_2`, respectively.
#' - Harmonizes category colors across both axes using a single legend.
#'
#' If `annotation_colors` is not supplied, a default palette is used. The resulting plot includes the heatmap, the sidebars, and a unified legend for variable categories.
#'
#' @return A `ggplot2` object assembled using `patchwork`, showing the correlation heatmap and annotation bars.
#'
#' @examples
#' \dontrun{
#' plot_heatmap_correlate_exposure(expom)
#' plot_heatmap_correlate_exposure(expom, exposure_cols = c("pm25", "no2", "age"))
#' plot_heatmap_correlate_exposure(expom, corr_threshold = 0.2)
#' }
#'
#' @export

plot_heatmap_correlate_exposure <- function(
    expomicset,
    exposure_cols = NULL,
    corr_threshold = NULL,
    annotation_colors = NULL,
    low = "#006666",
    mid = "white",
    high = "#8E0152",
    midpoint = 0) {
  require(ggplot2)
  require(patchwork)

  # Check that exposure correlation has been run
  if(!("exposure_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))) ) {
    stop("Exposure correlation has not been run. Please run `correlate_exposures()` first.")
  }

  correlation_df <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("exposure_correlation") |>
    purrr::pluck("correlation_table")

  if(!is.null(exposure_cols)){
    # Filter correlation_df based on exposure_cols
    correlation_df <- correlation_df |>
      dplyr::filter(var1 %in% exposure_cols & var2 %in% exposure_cols)
  }

  if(!is.null(corr_threshold)) {
    # Filter based on correlation threshold
    correlation_df <- correlation_df |>
      dplyr::filter(abs(correlation) >= corr_threshold)

    # Remove unnecessary rows if no correlations remain
    remaining_vars <- unique(c(correlation_df$var1, correlation_df$var2))

    correlation_df <- correlation_df |>
      dplyr::filter(var1 %in% remaining_vars & var2 %in% remaining_vars)
  }

  # Build correlation matrix
  cor_mat <- correlation_df |>
    dplyr::select(var1, var2, correlation) |>
    tidyr::pivot_wider(names_from = var2,
                       values_from = correlation) |>
    dplyr::mutate(across(everything(), ~replace_na(.x, 0))) |>
    tibble::column_to_rownames("var1") |>
    as.matrix()

  # Get variable order
  vars <- rownames(cor_mat)
  idx <- setNames(seq_along(vars), vars)

  # Melt correlation matrix and join with metadata
  df <- as.data.frame(as.table(cor_mat)) |>
    dplyr::rename(var1 = Var1,
                  var2 = Var2,
                  correlation = Freq) |>
    dplyr::left_join(correlation_df |>
                       dplyr::select(var1,
                                     var2,
                                     category_1,
                                     category_2),
              by = c("var1", "var2")) |>
    dplyr::mutate(
      var1_idx = idx[as.character(var1)],
      var2_idx = idx[as.character(var2)]
    ) |>
    dplyr::filter(var1_idx <= var2_idx) |>  # Keep upper-left triangle
    dplyr::mutate(
      var1 = factor(var1, levels = vars),
      var2 = factor(var2, levels = rev(vars))  # reverse for y-axis
    )

  # Unified color palette
  all_categories <- unique(c(correlation_df$category_1,
                             correlation_df$category_2)) |>
    as.character() |>
    unique() |>
    na.omit()

  # Check if annotation_colors is provided
  if(!is.null(annotation_colors)){
    cat_colors <- annotation_colors

    # Check if the length of annotation_colors matches the number of unique categories
    if(length(cat_colors) != length(all_categories)) {
      stop("Length of annotation_colors must match number of unique categories in the data.")
    }

  } else{
    cat_colors <- setNames(tidy_exp_pal[1:length(all_categories)], all_categories)
  }

  # Heatmap
  heatmap <- df |>
    ggplot(
      aes(x = var1,
          y = var2,
          fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  # Top annotation bar (x-axis)
  top_annot <- correlation_df |>
    dplyr::distinct(var1, category_1)

  top_bar <- top_annot |>
    ggplot(aes(
      x = var1,
      y = 1,
      fill = category_1)) +
    geom_tile(height = 0.01) +
    scale_fill_manual(name = "Category",
                      values = cat_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  # Left annotation bar (y-axis)
  left_annot <- correlation_df |>
    dplyr::distinct(var2, category_2)

  left_bar <- left_annot |>
    ggplot(aes(
      x = 1,
      y = var2,
      fill = category_2)) +
    geom_tile(width = 0.01) +
    scale_fill_manual(name = "Category",
                      values = cat_colors) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(hjust = 1),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )

  # Combine heatmap and annotations
  main_block <- heatmap / top_bar + plot_layout(heights = c(40, 2))
  left_stack <- left_bar / plot_layout(heights = c(40, 2))

  # Final layout
  (left_stack | main_block) +
    plot_layout(widths = c(2, 40), guides = "collect") &
    theme(legend.position = "right")
}
