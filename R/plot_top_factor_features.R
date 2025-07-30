#' Plot Top Features by Factor from Integration Results
#'
#' Visualizes the top loading features for each factor from multi-omics
#' integration results (e.g., MOFA, MCIA, DIABLO, RGCCA).
#'
#' @param expomicset A `MultiAssayExperiment` object containing
#' integration results in the `metadata` slot (must include `integration_results`).
#' @param factors Character vector of factors to include
#' (e.g., "Factor1", "Factor2"). If `NULL`, all factors are plotted.
#' @param top_n Integer specifying the number of top features to show
#'  per factor. Default is `5`.
#' @param facet_cols Optional color palette for facet strip backgrounds
#' (one per `exp_name`), used to distinguish factors.
#' @param exp_name_cols Optional color palette for experiment labels
#' in the plot (`exp_name`), passed to `scale_color_manual()`.
#' @param alpha Numeric value between 0 and 1 controlling the
#' transparency of facet strip background fill. Default is `0.5`.
#'
#' @details
#' This function supports the following integration methods:
#' - `"MOFA"`: Uses feature weights from MOFA2 (`get_weights()`).
#' - `"MCIA"`: Uses block loadings from MCIA (`@block_loadings`).
#' - `"DIABLO"`: Extracts block-specific loadings from `loadings`.
#' - `"RGCCA"`: Extracts block-specific loadings from `a`.
#'
#' For each factor, it:
#' - Selects the top `top_n` features by **absolute loading**.
#' - Creates a point-range plot showing the loading magnitude.
#' - Facets each factor with a customizable strip background.
#'
#' If palettes are not provided, defaults are chosen using
#' `ggpubr::get_palette()`.
#'
#' @return A `ggplot2` object with one facet per factor, showing the
#' top features and their loadings by experiment.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'    n_samples = 20,
#'    return_mae=TRUE
#'   )
#'
#' mae <- run_multiomics_integration(
#'       mae,
#'       method = "MCIA",
#'       n_factors = 3)
#'
#'
#' top_feature_p <- mae |>
#'   plot_top_factor_features()
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom MOFA2 get_weights
#' @importFrom purrr map2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate filter group_by arrange desc slice_head
#' bind_rows inner_join sym
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_bw
#' theme labs scale_color_manual element_text
#' @importFrom ggh4x facet_grid2 strip_themed elem_list_rect
#' @importFrom scales alpha
#' @export
plot_top_factor_features <- function(
    expomicset,
    feature_col="feature",
    factors=NULL,
    top_n=5,
    facet_cols=NULL,
    exp_name_cols=NULL,
    alpha=0.5){

  # require(ggplot2)

  method <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method
  result <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result

  loadings_df <- switch(
    method,

    "MOFA" = {
      MOFA2::get_weights(result) |>
        purrr::map2(
          names(MOFA2::get_weights(result)),
          ~as.data.frame(.x) |>
            tibble::rownames_to_column("feature") |>
            tidyr::pivot_longer(-feature,
                                names_to="factor",
                                values_to="loading") |>
            dplyr::mutate(abs_loading = abs(loading),
                          exp_name = .y)
        ) |>
        dplyr::bind_rows()
    },

    "MCIA" = {
      result@block_loadings |>
        purrr::map2(
          names(result@block_loadings),
          ~as.data.frame(.x) |>
            tibble::rownames_to_column("feature") |>
            tidyr::pivot_longer(-feature,
                                names_to="factor",
                                values_to="loading") |>
            dplyr::mutate(abs_loading = abs(loading),
                          exp_name = .y)
        ) |>
        dplyr::bind_rows()
    },

    "DIABLO" = {
      result$loadings |>
        purrr::map2(
          names(result$loadings),
          ~as.data.frame(.x) |>
            tibble::rownames_to_column("feature") |>
            tidyr::pivot_longer(-feature,
                                names_to="factor",
                                values_to="loading") |>
            dplyr::mutate(abs_loading = abs(loading),
                          exp_name = .y)
        ) |>
        dplyr::bind_rows() |>
        dplyr::mutate(factor = paste(exp_name,factor,sep = " "))
    },

    "RGCCA" = {
      result$a |>
        purrr::map2(
          names(result$a),
          ~as.data.frame(.x) |>
            tibble::rownames_to_column("feature") |>
            tidyr::pivot_longer(-feature,
                                names_to="factor",
                                values_to="loading") |>
            dplyr::mutate(abs_loading = abs(loading),
                          exp_name = .y)
        ) |>
        dplyr::bind_rows() |>
        dplyr::mutate(factor = paste(exp_name,factor,sep = " "))
    },

    stop("Method not supported.")
  )

  # Filter by factors if provided
  if (!is.null(factors)) {
    loadings_df <- loadings_df |>
      dplyr::filter(factor %in% factors)
  }

  # Select top features per factor
  df <- loadings_df |>
    dplyr::group_by(factor) |>
    dplyr::arrange(dplyr::desc(abs_loading)) |>
    dplyr::slice_head(n=top_n)

  # Map to the feature data
  df <- df |>
    inner_join(pivot_feature(expomicset),
               by=c("feature"=".feature",
                    "exp_name"=".exp_name"))

  # If no facet_cols provided, use default
  facet_cols <- facet_cols %||% tidy_exp_pal[
    seq_len(length(unique(df$factor)))
  ]

  # If no exp_name_cols provided, use default
  exp_name_cols <- exp_name_cols %||% rev(tidy_exp_pal)[
    seq_len(length(unique(df$exp_name)))
  ]


  # Create a plot of top features per factor
  features_per_factor_plot <- df |>
    ggplot(aes(
      x=abs_loading,
      y=reorder(!!dplyr::sym(feature_col), abs_loading),
      color=exp_name
    )) +
    geom_point(shape=18, size=5, alpha=0.5) +
    geom_segment(aes(
      x=0,
      xend=abs_loading,
      y=!!dplyr::sym(feature_col),
      yend=!!dplyr::sym(feature_col)),
      color="grey55") +
    theme_bw() +
    theme(
      strip.text.y = element_text(face="bold.italic", angle = 0),
      axis.text.y = element_text(face="italic")
    ) +
    ggh4x::facet_grid2(
      factor~.,
      scales = "free_y",
      space = "free_y",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = scales::alpha(facet_cols, alpha)
        )
      )
    ) +
    scale_color_manual(values=exp_name_cols) +
    labs(
      x="Absolute loading",
      y="",
      color="Experiment"
    )

  return(features_per_factor_plot)
}
