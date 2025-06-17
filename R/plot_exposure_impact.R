#' Visualize Exposure Impact Metrics from DEG Centrality Analysis
#'
#' This function generates a two-panel visualization of the results from
#' \code{\link{run_exposure_impact}}. It displays (1) the degree centrality
#' of features associated with each exposure as a boxplot, and (2) a stacked
#' bar plot summarizing how many differentially expressed features are associated
#' with exposures.
#'
#' @param expomicset A \code{MultiAssayExperiment} object containing the output
#'   from \code{\link{run_exposure_impact}} in its metadata.
#' @param min_per_group Integer. Minimum number of features required for an
#'   exposure group to be included in the boxplot. Default is \code{5}.
#' @param facet_cols Optional vector of colors for the exposure category facets.
#'   If \code{NULL}, defaults to \code{ggsci::pal_jco()} colors.
#' @param bar_cols Optional vector of colors for the bar plot groups. If \code{NULL},
#'   defaults to \code{ggsci::pal_uchicago()} colors.
#' @param alpha Numeric between 0 and 1. Transparency of facet background fill.
#'   Default is \code{0.5}.
#' @param ncol Integer. Number of columns in the combined patchwork layout.
#'   Default is \code{2}.
#' @param nrow Integer. Number of rows in the combined patchwork layout.
#'   Default is \code{1}.
#' @param heights Numeric vector. Relative heights of the patchwork panels.
#'   Default is \code{c(1, 1)}.
#' @param widths Numeric vector. Relative widths of the patchwork panels.
#'   Default is \code{c(2, 1)}.
#'
#' @return A \code{ggplot} object composed with \code{patchwork} that includes:
#' \enumerate{
#'   \item A boxplot of degree centrality values for features associated with each exposure.
#'   \item A barplot showing the number of DEGs associated and not associated with exposures.
#' }
#'
#' @details
#' This function assumes that \code{\link{run_exposure_impact}} has already been run and that
#' the resulting exposure impact metrics are stored in the \code{metadata(expomicset)} list
#' under the \code{"exposure_impact"} entry.
#'
#' The left panel (boxplot) stratifies centrality by exposure and facet category.
#' The right panel (barplot) summarizes the total and associated DEG counts.
#'
#' @seealso \code{\link{run_exposure_impact}}, \code{\link[patchwork]{wrap_plots}}
#'
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

  exposure_impact_deg_n <-
    expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("network") |>
    purrr::pluck("exposure_impact") |>
    purrr::pluck(feature_type) |>
    purrr::pluck("deg_association_summary")

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
        exposure_impact_deg_n |>
          dplyr::pull(exp_name)
      )))
  }

  # Create a boxplot of the degree centrality of genes that are associated
  # with a particular exposure
  boxplot <- exposure_impact_degree |>
    dplyr::filter(n_features >min_per_group) |>
    ggplot(aes(
      x=reorder(exposure,-mean_centrality),
      y=centrality,
      fill=category,
      color=category
    ))+
    geom_boxplot(alpha=0.7) +
    geom_jitter(alpha=0.1)+
    ggpubr::rotate_x_text(angle=45)+
    theme_minimal()+
    #ggpubr::theme_pubclean()+
    scale_fill_manual(values = facet_cols) +
    ggpubr::rotate_x_text(angle=45)+
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
    scale_color_manual(values=facet_cols)+
    theme(strip.text.x = element_text(face = "bold.italic"),
          legend.position = "none")+
    labs(
      x="",
      y="Degree",
    )

  barplot <- exposure_impact_deg_n |>
    dplyr::select(-percent_assoc_exp) |>
    dplyr::mutate(left_over= total_deg-n_assoc_exp) |>
    dplyr::select(-total_deg) |>
    tidyr::pivot_longer(-exp_name,
                        names_to = "group",
                        values_to = "value") |>
    dplyr::mutate(group=dplyr::case_when(
      group == "n_assoc_exp" ~ "Associated with Exposures",
      group == "left_over" ~ "Not Associated with Exposures"
    )) |>
    ggplot(aes(
      x=value,
      y=reorder(exp_name,value),
      fill=group
    ))+
    geom_bar(stat="identity")+
    ggpubr::rotate_x_text(angle=45)+
    ggpubr::theme_pubr(legend="right")+
    scale_fill_manual(values = bar_cols)+
    labs(
      x="Number of Features",
      y="",
      fill=""
    )

  final_plot <- patchwork::wrap_plots(
    boxplot,
    barplot,
    ncol = ncol,
    nrow = nrow)+
    patchwork::plot_layout(
      heights = heights,
      widths = widths
    )


  return(final_plot)
}
