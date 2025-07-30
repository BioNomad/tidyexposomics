#' Plot Exposure Distributions by Category or Group
#'
#' Visualizes exposure variable distributions using **boxplots**
#' or **ridge plots**.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure data.
#' @param exposure_cat A character string or vector specifying exposure
#' category names (from `codebook$category`) to include.
#'  Use `"all"` to include all exposures.
#' @param exposure_cols Optional character vector specifying exact
#' exposure variables to plot.
#' @param group_by A string specifying the column in `colData(expomicset)`
#' used to fill the plot (e.g., `"sex"`). Defaults to `NULL`, in which case
#' exposures are colored by `category`.
#' @param plot_type Type of plot: `"boxplot"` (default) or `"ridge"`.
#' @param alpha Transparency level for background facet color strips.
#' Default is `0.5`.
#' @param panel_sizes A numeric vector passed to `ggh4x::force_panelsizes()`
#' for controlling facet widths or heights.
#' @param title Plot title. Default is `"Exposure Levels by Category"`.
#' @param xlab X-axis label. Default is an empty string.
#' @param ylab Y-axis label. Default is an empty string.
#' @param facet_cols Optional vector of colors to use as background for
#'  facet categories. If `NULL`, a default palette is used.
#' @param group_cols Optional named vector of colors for `group_by` levels.
#' If `NULL`, a default palette is used.
#' @param fill_lab Legend title for the fill aesthetic
#' (e.g., `"Sex"` or `"Exposure Group"`). Default is `""`.
#'
#' @details
#' This function:
#' - Filters exposure data based on category or selected columns.
#' - Merges variable metadata from `metadata(expomicset)$codebook`.
#' - Supports either **boxplot** (vertical distributions per variable)
#' or **ridgeplot** (horizontal density plots per variable).
#' - If `group_by` is specified, that variable defines the plot fill color;
#' otherwise, the fill is based on exposure `category`.
#' - Facets by `category` using `ggh4x::facet_grid2()`
#' with color-coded strip backgrounds.
#'
#' @return A `ggplot2` object showing exposure distributions,
#' optionally grouped.
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 10,
#'   return_mae=TRUE
#'   )
#'
#' # plot exposure data
#' exposure_plot <- mae |>
#'   plot_exposures(
#'     exposure_cols = c("exposure_pm25","exposure_no2")
#'   )
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter
#' geom_tile scale_color_manual scale_fill_manual scale_y_log10
#' scale_x_log10 labs theme_minimal theme element_text guides guide_legend
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter select pull inner_join all_of where distinct mutate
#' @importFrom purrr pluck
#' @importFrom ggh4x facet_grid2 strip_themed elem_list_rect force_panelsizes
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggpubr theme_pubclean
#'
#' @export
plot_exposures <- function(
    expomicset,
    exposure_cat = "all",
    exposure_cols = NULL,
    group_by = NULL,
    plot_type = "boxplot",
    alpha=0.3,
    panel_sizes=rep(1,100),
    title = "Exposure Levels by Category",
    xlab = "",
    ylab = "",
    facet_cols = NULL,
    group_cols = NULL,
    box_width = 0.1,
    fill_lab = ""
){
  # require(ggplot2)

  # Extract exposure data
  exposure_data <- expomicset |>
    pivot_sample()

  # Extract variable description file
  des <- MultiAssayExperiment::metadata(expomicset) |>
    purrr::pluck("codebook")

  # Filter by exposure category if specified
  if (exposure_cat == "all") {
    exposure_data <- exposure_data
  } else {
    vars_to_keep <- des |>
      dplyr::filter(category %in% exposure_cat) |>
      dplyr::pull(variable)
    cols_to_keep <- c(".sample", vars_to_keep)
    if (!is.null(group_by)) cols_to_keep <- c(cols_to_keep, group_by)
    cols_to_keep <- intersect(cols_to_keep, colnames(exposure_data))
    exposure_data <- exposure_data[, cols_to_keep]
  }

  # If specific columns are provided, filter to those
  if (!is.null(exposure_cols)) {
    exposure_data <- exposure_data[, c(".sample",group_by, exposure_cols)]
  }

  # Ensure only numeric exposure variables + .sample are retained
  numeric_cols <- exposure_data |>
    dplyr::select(dplyr::where(is.numeric)) |>
    colnames()

  exposure_data <- exposure_data |>
    dplyr::select(dplyr::all_of(c(".sample", group_by,numeric_cols)))


  # Pivot to long format and join with variable metadata
  sample_metadata <- exposure_data |>
    tidyr::pivot_longer(
      cols = -c(.sample, group_by),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::inner_join(des, by = "variable")

  # set facet colors
  if(!is.null(facet_cols)){
    facet_cols <- facet_cols
  } else{
    # facet_cols <- tidy_exp_pal[1:length(unique(sample_metadata$category))]
    facet_cols <- tidy_exp_pal[
      seq_len(length(unique(sample_metadata$category)))
    ]
  }

  if(!is.null(group_by)){
    group_var <- group_by
  } else {
    group_var <- "category"
  }

  if(!is.null(group_cols)){
    group_cols <- group_cols
  } else{
    # group_cols <- tidy_exp_pal[1:length(unique(sample_metadata[[group_var]]))]
    group_cols <- tidy_exp_pal[
      seq_len(length(unique(sample_metadata[[group_var]])))
    ]
  }

  if( plot_type == "boxplot"){

    legend_logic <- ifelse(!is.null(group_by),"right","none")

    # Create boxplot
    sample_metadata |>
      dplyr::filter(value>0) |>
      ggplot(aes(
        x = variable,
        y = value,
        fill = !!sym(group_var),
        color = !!sym(group_var))) +
      #geom_violin()+
      geom_boxplot(alpha=0.4) +
      geom_jitter(alpha=0.03)+
      ggh4x::facet_grid2(
        ~category,
        scales = "free_x",
        space = "free_x",
        strip = ggh4x::strip_themed(
          background_x = ggh4x::elem_list_rect(
            fill = scales::alpha(
              facet_cols,
              alpha
            )
          )
        )
      ) +
      scale_color_manual(values = group_cols) +
      scale_fill_manual(values = group_cols) +
      ggh4x::force_panelsizes(cols = panel_sizes)+
      labs(
        title = title,
        x = xlab,
        y = ylab,
        fill = fill_lab
      ) +
      #ggpubr::theme_pubclean() +
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 55, hjust = 1),
            legend.position = legend_logic,
            strip.text.x = element_text(face = "bold.italic"),
            plot.title = element_text(face = "bold.italic"))+
      guides(color = "none")+
      scale_y_log10()

  } else if (plot_type == "ridge"){
    # Create boxplot
    sample_metadata |>
      filter(value>0) |>
      ggplot(aes(
        x = value,
        y = variable,
        fill = !!sym(group_var))) +
      ggridges::geom_density_ridges() +
      ggh4x::facet_grid2(
        category~.,
        scales = "free_y",
        space = "free_y",
        strip = ggh4x::strip_themed(
          background_y = ggh4x::elem_list_rect(
            fill = scales::alpha(
              facet_cols,
              alpha
            )
          )
        )
      ) +
      scale_fill_manual(values = group_cols)+
      ggh4x::force_panelsizes(cols = panel_sizes)+
      labs(
        title = title,
        x = xlab,
        y = ylab,
        fill = fill_lab
      ) +
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = ifelse(!is.null(group_var),"right","none"),
            strip.text.y = element_text(angle=0,face = "bold.italic"),
            plot.title = element_text(face = "bold.italic"))+
      scale_x_log10()
  } else {
    stop("Unsupported plot type.")
  }


}


