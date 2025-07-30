#' Plot Summary of Factor Contributions from Multi-Omics Integration
#'
#' Generates a summary plot of factor contributions from multi-omics
#' integration results stored in a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing integration
#' results in
#'   `metadata(expomicset)$multiomics_integration$integration_results`.
#' @param low Color for low values in the fill gradient.
#' Default is `"#006666"`.
#' @param mid Color for midpoint in the fill gradient.
#' Default is `"white"`.
#' @param high Color for high values in the fill gradient.
#' Default is `"#8E0152"`.
#' @param midpoint Midpoint value for the gradient color scale.
#' Default is `0.5`.
#'
#' @details
#' This function visualizes factor contributions based on the integration method:
#'
#' - **MOFA**: Variance explained per factor and view.
#' - **MCIA**: Block score weights per omic.
#' - **DIABLO**: Mean absolute sample score per omic and factor
#' (from block-specific variates).
#' - **RGCCA**: Mean absolute sample score per omic and factor
#' (from aligned block scores).
#'
#' The color gradient can be customized using the `low`, `mid`, `high`,
#' and `midpoint` parameters.
#'
#' @return A `ggplot` object showing factor contributions based on
#'  the integration method.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'   n_samples = 20,
#'   return_mae=TRUE
#'   )
#'
#' mae <- run_multiomics_integration(
#'       mae,
#'       method = "MCIA",
#'       n_factors = 3)
#'
#' factor_sum_plot <- mae |>
#'   plot_factor_summary()
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs
#'  theme_minimal element_text
#' @importFrom ggpubr theme_pubr rotate_x_text
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate bind_rows select distinct
#' @importFrom purrr discard imap map
#' @export
plot_factor_summary <- function(
    expomicset,
    low = "#006666",
    mid = "white",
    high = "#8E0152",
    midpoint = 0.5
) {
  # require("ggplot2")
  # require("ggpubr")

  integration <- MultiAssayExperiment::metadata(expomicset) |>
    purrr::pluck("multiomics_integration",
                 "integration_results")

  if (is.null(integration) || is.null(integration$result)) {
    stop("Integration results not found.
         Please run `run_multiomics_integration()` first.")
  }

  method <- integration$method
  result <- integration$result

  factor_contrib_plot <- switch(
    method,

    "MOFA" = {
      MOFA2::plot_variance_explained(
        result,
        x = "view",
        y = "factor"
      ) +
        ggpubr::rotate_x_text(45) +
        scale_fill_gradient2(low = low,
                             mid = mid,
                             high = high,
                             midpoint = midpoint)
    },

    "MCIA" = {
      result@block_score_weights |>
        as.data.frame() |>
        tibble::rownames_to_column("omic") |>
        tidyr::pivot_longer(!omic,
                            names_to = "factor",
                            values_to = "weight") |>
        dplyr::mutate(
          factor = gsub("V", "", factor),
          factor = factor(as.numeric(factor),
                          levels = sort(unique(as.numeric(factor))))
        ) |>
        ggplot(aes(x = factor,
                   y = omic,
                   fill = weight)) +
        geom_tile() +
        ggpubr::theme_pubr(legend = "right") +
        scale_fill_gradient2(low = low,
                             mid = mid,
                             high = high,
                             midpoint = midpoint) +
        labs(x = "Factor",
             y = NULL,
             fill = "Weight")
    },

    # "MCCA" = {
    #   result$sample_scores |>
    #     purrr::map(~ apply(.x, 2, function(x) mean(abs(x))) |>
    #                  as.data.frame() |>
    #                  tibble::rownames_to_column("factor") |>
    #                  setNames(c("factor", "weight"))) |>
    #     dplyr::bind_rows(.id = "omic") |>
    #     dplyr::mutate(
    #       factor = factor(as.numeric(factor),
    #                       levels = sort(unique(as.numeric(factor))))
    #     ) |>
    #     ggplot(aes(x = factor, y = omic, fill = weight)) +
    #     geom_tile() +
    #     ggpubr::theme_pubr(legend = "right") +
    #     scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = midpoint) +
    #     labs(x = "Factor", y = "", fill = "Avg |Score|")
    # },

    "DIABLO" = {
      result$variates |>
        (\(lst) lst[names(lst) != "Y"])() |>
        purrr::discard(~ is.null(.x) || is.character(.x) || is.factor(.x)) |>
        purrr::imap(~ apply(.x, 2, function(x) mean(abs(x))) |>
                      as.data.frame() |>
                      tibble::rownames_to_column("factor") |>
                      setNames(c("factor", "weight"))) |>
        dplyr::bind_rows(.id = "omic") |>
        dplyr::mutate(
          factor = gsub("comp", "", factor),
          factor = factor(as.numeric(factor),
                          levels = sort(unique(as.numeric(factor))))
        ) |>
        ggplot(aes(x = factor, y = omic, fill = weight)) +
        geom_tile() +
        ggpubr::theme_pubr(legend = "right") +
        scale_fill_gradient2(low = low,
                             mid = mid,
                             high = high,
                             midpoint = midpoint) +
        labs(x = "Factor",
             y = "",
             fill = "Avg. |Score|")
    },

    "RGCCA" = {
      result$Y |>
        purrr::imap(~ apply(.x, 2, function(x) mean(abs(x))) |>
                      as.data.frame() |>
                      tibble::rownames_to_column("factor") |>
                      setNames(c("factor", "weight"))) |>
        dplyr::bind_rows(.id = "omic") |>
        dplyr::mutate(
          factor = gsub("comp", "", factor),
          factor = factor(as.numeric(factor),
                          levels = sort(unique(as.numeric(factor))))
        ) |>
        ggplot(aes(x = factor, y = omic, fill = weight)) +
        geom_tile() +
        ggpubr::theme_pubr(legend = "right") +
        scale_fill_gradient2(low = low,
                             mid = mid,
                             high = high,
                             midpoint = midpoint) +
        labs(x = "Factor",
             y = "",
             fill = "Avg. |Score|")
    },

    stop("Integration method not supported in plot_factor_summary().")
  )

  return(factor_contrib_plot)
}
