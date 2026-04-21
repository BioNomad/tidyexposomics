#' Plot Exposure-Omics Associations
#'
#' @description
#' Plots the number of significant exposure-omics associations, grouped
#' either by exposure or the exposure category.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing association results.
#' @param plot_type Character. One of `"exposures"` or `"category"`. Controls
#'   whether associations are summarized per exposure or per exposure category.
#'   Defaults to `"exposures"`.
#' @param pval_col Character. Name of the column used for p-value filtering.
#'   Defaults to `"p_adjust"`.
#' @param pval_thresh Numeric. Significance threshold applied to `pval_col`.
#'   Rows with values below this threshold are retained. Defaults to `0.05`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # run exposure-omics association
#' mae <- mae |>
#'     run_exposure_omics_association(
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         covariates = c("age", "sex")
#'     )
#'
#' plot_exposure_omics_association(
#'     exposomicset = mae,
#'     plot_type = "exposures"
#' )
#'
#' @export
plot_exposure_omics_association <- function(
  exposomicset,
  plot_type = c("exposures", "category"),
  pval_col = "p_adjust",
  pval_thresh = 0.05
) {
    plot <- switch(plot_type,
        "exposures" = .plot_exp_omic_exposure(
            exposomicset,
            pval_col = pval_col,
            pval_thresh = pval_thresh
        ),
        "category" = .plot_exp_omic_category(
            exposomicset,
            pval_col = pval_col,
            pval_thresh = pval_thresh
        )
    )
    return(plot)
}

#' Plot Significant Omics Associations per Exposure
#'
#' @description
#' Internal helper that counts and plots the number of significant omics
#' associations for each individual exposure.
#'
#' @inheritParams plot_exposure_omics_association
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
.plot_exp_omic_exposure <- function(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05
) {
    extract_results(
        exposomicset = exposomicset,
        result = "association"
    ) |>
        purrr::pluck(
            "exposure_omics",
            "results_df"
        ) |>
        filter(!!sym(pval_col) < pval_thresh) |>
        group_by(exposure) |>
        dplyr::count() |>
        ggplot(aes(
            x = n,
            y = reorder(exposure, n),
            fill = n
        )) +
        geom_col(linewidth = 0.005, color = "black") +
        scale_fill_gradient(
            low = "#fdeff9",
            high = "#03001e"
        ) +
        ggpubr::theme_pubr(legend = "right") +
        labs(
            x = "No. of Omics Associations",
            y = "",
            fill = "No. of Omics Associations",
            title = "No. of Omics Associations per Exposure"
        ) +
        theme(
            plot.title = element_text(face = "bold.italic")
        )
}

#' Plot Significant Omics Associations per Exposure Category
#'
#' @description
#' Internal helper that counts and plots the number of significant omics
#' associations for each exposure category.
#'
#' @inheritParams plot_exposure_omics_association
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
.plot_exp_omic_category <- function(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05
) {
    extract_results(
        exposomicset = exposomicset,
        result = "association"
    ) |>
        purrr::pluck(
            "exposure_omics",
            "results_df"
        ) |>
        filter(!!sym(pval_col) < pval_thresh) |>
        group_by(category) |>
        dplyr::count() |>
        ggplot(aes(
            x = n,
            y = reorder(category, n),
            fill = n
        )) +
        geom_col(linewidth = 0.005, color = "black") +
        scale_fill_gradient(
            low = "#fdeff9",
            high = "#03001e"
        ) +
        ggpubr::theme_pubr(legend = "right") +
        labs(
            x = "No. of Omics Associations",
            y = "",
            fill = "No. of Omics Associations",
            title = "No. of Omics Associations per Category"
        ) +
        theme(
            plot.title = element_text(face = "bold.italic")
        )
}
