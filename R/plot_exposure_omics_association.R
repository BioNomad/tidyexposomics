#' Plot Exposure-Omics Associations
#'
#' @description
#' Plots the number of significant exposure-omics associations, grouped
#' by exposure, exposure category, or omics experiment.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing association results.
#' @param plot_type Character. One of `"exposures"`, `"category"`, or
#'   `"exp_name"`. Controls whether associations are summarized per exposure,
#'   per exposure category, or per omics experiment (assay). Defaults to
#'   `"exposures"`.
#' @param pval_col Character. Name of the column used for p-value filtering.
#'   Defaults to `"p_adjust"`.
#' @param pval_thresh Numeric. Significance threshold applied to `pval_col`.
#'   Rows with values below this threshold are retained. Defaults to `0.05`.
#' @param estimate_thresh Numeric. Minimum absolute effect size. Rows with
#'   `abs(estimate)` greater than or equal to this value are retained.
#'   Defaults to `0` (no filtering on estimate).
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
#'         formula = ~ .feature + age + sex,
#'         exposures = c("exposure_pm25", "exposure_no2")
#'     )
#'
#' plot_exposure_omics_association(
#'     exposomicset = mae,
#'     plot_type = "exp_name",
#'     estimate_thresh = 0.1
#' )
#'
#' @export
plot_exposure_omics_association <- function(
  exposomicset,
  plot_type = c("exposures", "category", "exp_name"),
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  estimate_thresh = 0
) {
    plot_type <- match.arg(plot_type)

    switch(plot_type,
        "exposures" = .plot_exp_omic_exposure(
            exposomicset,
            pval_col = pval_col,
            pval_thresh = pval_thresh,
            estimate_thresh = estimate_thresh
        ),
        "category" = .plot_exp_omic_category(
            exposomicset,
            pval_col = pval_col,
            pval_thresh = pval_thresh,
            estimate_thresh = estimate_thresh
        ),
        "exp_name" = .plot_exp_omic_exp_name(
            exposomicset,
            pval_col = pval_col,
            pval_thresh = pval_thresh,
            estimate_thresh = estimate_thresh
        )
    )
}

#' Get Filtered Exposure-Omics Association Results
#'
#' @description
#' Internal helper that extracts the exposure-omics results and keeps rows
#' passing both the p-value and effect size thresholds.
#'
#' @inheritParams plot_exposure_omics_association
#'
#' @return A tibble of filtered association results.
#'
#' @keywords internal
.get_exp_omic_results <- function(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  estimate_thresh = 0
) {
    extract_results(
        exposomicset = exposomicset,
        result = "association"
    ) |>
        purrr::pluck(
            "exposure_omics",
            "results_df"
        ) |>
        filter(
            !!sym(pval_col) < pval_thresh,
            abs(estimate) >= estimate_thresh
        )
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
  pval_thresh = 0.05,
  estimate_thresh = 0
) {
    .get_exp_omic_results(
        exposomicset,
        pval_col = pval_col,
        pval_thresh = pval_thresh,
        estimate_thresh = estimate_thresh
    ) |>
        dplyr::count(exposure) |>
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
  pval_thresh = 0.05,
  estimate_thresh = 0
) {
    .get_exp_omic_results(
        exposomicset,
        pval_col = pval_col,
        pval_thresh = pval_thresh,
        estimate_thresh = estimate_thresh
    ) |>
        dplyr::count(category) |>
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

#' Plot Significant Associations per Omics Experiment
#'
#' @description
#' Internal helper that counts and plots the number of significant
#' exposure-omics associations for each omics experiment (assay).
#'
#' @inheritParams plot_exposure_omics_association
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
.plot_exp_omic_exp_name <- function(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  estimate_thresh = 0
) {
    .get_exp_omic_results(
        exposomicset,
        pval_col = pval_col,
        pval_thresh = pval_thresh,
        estimate_thresh = estimate_thresh
    ) |>
        dplyr::count(exp_name) |>
        ggplot(aes(
            x = n,
            y = reorder(exp_name, n),
            fill = n
        )) +
        geom_col(linewidth = 0.005, color = "black") +
        scale_fill_gradient(
            low = "#fdeff9",
            high = "#03001e"
        ) +
        ggpubr::theme_pubr(legend = "right") +
        labs(
            x = "No. of Associations",
            y = "",
            fill = "No. of Associations",
            title = "No. of Associations per Omics Experiment"
        ) +
        theme(
            plot.title = element_text(face = "bold.italic")
        )
}
