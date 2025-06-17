#' Plot Association Results (Unified Forest Plot)
#'
#' Generates a forest plot for association results from any source stored in the metadata of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object with association results in metadata.
#' @param source One of `"omics"`, `"exposures"`, `"factors"`, `"go_pcs"`.
#' @param filter_col Column to use for significance filtering. Default is `"p.value"`.
#' @param filter_thresh Numeric threshold for `filter_col`. Default is `0.05`.
#' @param direction_filter One of `"all"`, `"up"`, `"down"`. Default is `"all"`.
#' @param facet Logical; if `TRUE` and `source == "go_pcs"`, apply nested faceting. Default is `FALSE`.
#' @param nrow Facet row layout if faceting is enabled. Default is `1`.
#' @param subtitle Optional plot subtitle. If `NULL`, generated from covariates.
#'
#' @return A `ggplot2` forest plot of significant associations.
#'
#' @export
plot_association <- function(
    expomicset,
    source = c("omics", "exposures", "factors", "go_pcs"),
    terms = NULL,
    filter_col = "p.value",
    filter_thresh = 0.05,
    direction_filter = "all",
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
    stop("Association results for source '", source, "' not found in metadata. Run `run_association()` first.")
  }

  assoc_result <- assoc_metadata[[assoc_key]]
  results_df <- assoc_result$results_df
  covariates <- assoc_result$covariates

  # Filter significant associations
  df <- results_df |>
    dplyr::filter(!!rlang::sym(filter_col) < filter_thresh) |>
    dplyr::mutate(direction = dplyr::case_when(
      estimate > 0 ~ "up",
      estimate < 0 ~ "down",
      .default = "ns"
    ))

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

  if(!is.null(terms)){
    df <- dplyr::filter(df, term %in% terms)
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
      plot.subtitle = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::labs(
      x = "Effect size",
      y = "",
      title = paste0(tools::toTitleCase(source), "-Outcome Associations"),
      subtitle = subtitle
    )

  if (facet && source == "go_pcs") {
    if (!requireNamespace("ggh4x", quietly = TRUE)) {
      warning("Faceting requires the `ggh4x` package. Install it to enable nested facets.")
    } else {
      p <- p +
        ggh4x::facet_nested_wrap(facet_formula, nrow = nrow, scales = "free_x") +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic"))
    }
  }

  return(p)
}
