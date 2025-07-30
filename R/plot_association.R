#' Plot Association Results (Unified Forest Plot)
#'
#' Generates a forest plot for association results from any source stored
#' in the metadata of a `MultiAssayExperiment` object.
#' Supports faceting and visual augmentation with R^2 tiles when available.
#'
#' @param expomicset A `MultiAssayExperiment` object containing association
#' results in metadata.
#' @param source Character string indicating the association source.
#' One of `"omics"`, `"exposures"`, `"factors"`, or `"go_pcs"`.
#' @param terms Optional character vector of term names to subset the plot to.
#'  Default is `NULL` (include all).
#' @param filter_col Column used to assess statistical significance
#' (default: `"p.value"`).
#' @param filter_thresh Numeric threshold applied to `filter_col`
#' (default: `0.05`).
#' @param direction_filter Direction of associations to retain.
#' One of `"all"` (default), `"up"`, or `"down"`.
#' @param add_r2_tile Logical; if `TRUE`, includes a tile plot for `r2_col`
#'  (default: `TRUE`).
#' @param r2_col Column used for coloring the tile plot (default: `"adj_r2"`).
#' @param facet Logical; if `TRUE` and `source == "go_pcs"`,
#'  apply nested faceting by experiment and GO cluster (default: `FALSE`).
#' @param nrow Integer; number of rows for facet layout if enabled (default: `1`).
#' @param subtitle Optional subtitle for the plot. If `NULL`,
#' automatically generated from covariates used in the model.
#'
#' @return A `ggplot2` object: either a single forest plot or a composite plot
#'  with an R^2 tile strip.
#'
#' @details
#' This function visualizes effect size estimates and confidence intervals from
#' association analyses. It allows filtering by
#' direction (`"up"` for positive, `"down"` for negative) and by significance.
#'  For `source = "go_pcs"`, it supports special
#' formatting by splitting term labels into nested facets.
#'
#' The R^2 tile (if enabled) adds a side heatmap indicating model fit for
#' each association. This can be useful for model diagnostics.
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 10,
#'   return_mae=TRUE
#'   )
#'
#' # run association tests
#' mae <- mae |>
#'   run_association(
#'   source = "exposures",
#'   feature_set = c("exposure_pm25","exposure_no2"),
#'   outcome = "smoker",
#'   covariates = c("age"),
#'   family = "binomial"
#'   )
#'
#' assoc_plot <- mae |>
#'   plot_association(
#'     source = "exposures"
#'   )
#'
#' @export
plot_association <- function(
    expomicset,
    source = c("omics", "exposures", "factors", "go_pcs"),
    terms = NULL,
    filter_col = "p.value",
    filter_thresh = 0.05,
    direction_filter = "all",
    add_r2_tile = TRUE,
    r2_col = "adj_r2",
    facet = FALSE,
    nrow = 1,
    subtitle = NULL
) {
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  requireNamespace("tidyr")

  source <- match.arg(source)
  assoc_key <- paste0("assoc_", source)

  # Check metadata
  assoc_metadata <- MultiAssayExperiment::metadata(expomicset)$association
  if (is.null(assoc_metadata) || !assoc_key %in% names(assoc_metadata)) {
    stop("Association results for source '",
         source,
         "' not found in metadata. Run `run_association()` first.")
  }

  assoc_result <- assoc_metadata[[assoc_key]]
  results_df <- assoc_result$results_df
  covariates <- assoc_result$covariates

  if (!is.null(terms)) {
    df <- dplyr::filter(results_df, term %in% terms)
  } else {
    df <- results_df
  }

  # Mark significance
  df <- df |>
    dplyr::mutate(
      is_significant = !!rlang::sym(filter_col) < filter_thresh,
      direction = dplyr::case_when(
        (estimate > 0 & (is_significant == TRUE)) ~ "up",
        (estimate < 0 & (is_significant == TRUE)) ~ "down",
        .default = "ns"
      )
    )

  # Apply direction filter
  if (direction_filter == "up") {
    df <- dplyr::filter(df, direction == "up")
  } else if (direction_filter == "down") {
    df <- dplyr::filter(df, direction == "down")
  }

  # Extract y-axis variable and handle GO PC formatting
  if (source == "go_pcs") {
    df <- tidyr::separate(
      df,
      col = "term",
      into = c("PC", "exp_name", "Cluster", "go_group"),
      sep = "/",
      remove = FALSE
    )
    df$label <- df$go_group
    facet_formula <- exp_name + Cluster ~ .
  } else {
    df$label <- df$term
    facet_formula <- NULL
  }

  # Only proceed if add_r2_tile is TRUE and column exists
  if (add_r2_tile && r2_col %in% names(df)) {

    leg_label <- ifelse(r2_col == "r2",
                        expression("R"^2),
                        expression("Adj. R"^2))
    tile_plot <- df |>
      ggplot2::ggplot(ggplot2::aes(
                        x = 1,
                        y = reorder(label, estimate),
                        fill = !!rlang::sym(r2_col))) +
      ggplot2::geom_tile() +
      scale_fill_gradient(
        low = "#e8e1eb",
        high = "#3f1b4f"
      )+
      ggplot2::labs(x = NULL, y = NULL, fill = r2_col) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.text.y = element_blank()
      )+
      labs(fill=leg_label)
  }


  # Auto-subtitle
  if (is.null(subtitle)) {
    if (!is.null(covariates) && length(covariates) > 0) {
      subtitle <- paste("Covariates:", paste(covariates, collapse = ", "))
    }
  }

  # Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = estimate,
    y = reorder(label, estimate),
    color = direction
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = estimate - std.error, xmax = estimate + std.error),
      color = "grey55", height = 0.2
    ) +
    ggplot2::geom_point(shape = 18, size = 5, alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
      up = "#8E0152",
      down = "#006666",
      ns = "grey55"
    )) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.subtitle = ggplot2::element_text(face = "italic"),
      plot.title = ggplot2::element_text(face = "bold.italic"),
    ) +
    ggplot2::labs(
      x = "Effect size",
      y = "",
      title = paste0(tools::toTitleCase(source), "-Outcome Associations"),
      subtitle = subtitle
    )

  if (facet && source == "go_pcs") {
    if (!requireNamespace("ggh4x", quietly = TRUE)) {
      warning("Faceting requires the `ggh4x` package.
              Install it to enable nested facets.")
    } else {
      p <- p +
        ggh4x::facet_nested_wrap(
          facet_formula,
          nrow = nrow,
          scales = "free_x") +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"))
    }
  }

  if (add_r2_tile && exists("tile_plot")) {
    return(patchwork::wrap_plots(
      p + ggplot2::theme(axis.text.y = element_text()),
      tile_plot,
      widths = c(3, 0.5)))
  } else {
    return(p)
  }

  return(p)
}
