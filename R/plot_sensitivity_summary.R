#' Plot Sensitivity Analysis Summary
#'
#' Generates a ridge plot and bar chart summarizing feature stability
#' scores across assays.
#'
#' @param expomicset A `MultiAssayExperiment` object containing sensitivity
#'  analysis results
#' in `metadata(expomicset)$sensitivity_analysis`.
#' @param stability_score_thresh A numeric threshold for stability scores.
#' Default is `NULL`,
#' which uses the threshold stored in
#' `metadata(expomicset)$sensitivity_analysis$score_thresh`.
#' @param stability_metric A character string specifying which
#'  stability metric to plot (e.g., "stability_score",
#'   "logp_weighted_score"). Default is "stability_score".
#' @param title A character string specifying the title of the ridge plot.
#'  Default is
#' "Distribution of Stability Scores".
#'
#' @details
#' This function:
#' - Extracts feature stability scores from `metadata(expomicset)$sensitivity_analysis$feature_stability`.
#' - Displays a **ridge plot** of stability score distributions per assay.
#' - Displays a **bar chart** of the number of features per assay.
#' - Prints the number of features with stability scores above the threshold.
#'
#' @return A `patchwork` object combining a ridge plot and a bar chart.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#'
#' # Run differential abundance
#' mae <- run_differential_abundance(
#'     expomicset = mae,
#'     formula = ~ smoker + sex,
#'     abundance_col = "counts",
#'     method = "limma_voom",
#'     action = "add"
#' )
#'
#' # Run the sensitivity analysis
#' mae <- run_sensitivity_analysis(
#'     expomicset = mae,
#'     base_formula = ~ smoker + sex,
#'     methods = c("limma_voom"),
#'     scaling_methods = c("none"),
#'     covariates_to_remove = "sex",
#'     pval_col = "P.Value",
#'     logfc_col = "logFC",
#'     pval_threshold = 0.05,
#'     logFC_threshold = 0,
#'     bootstrap_n = 3,
#'     action = "add"
#' )
#'
#' # create the sensitivity summary plot
#' sens_sum_p <- mae |>
#'     plot_sensitivity_summary()
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom dplyr group_by summarise arrange left_join mutate n
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment scale_fill_manual
#'   scale_color_manual theme labs geom_vline element_text
#' @importFrom ggpubr theme_pubr
#' @export
plot_sensitivity_summary <- function(
    expomicset,
    stability_score_thresh = NULL,
    stability_metric = "stability_score",
    title = "Distribution of Stability Scores") {
    # require(ggplot2)
    # require(patchwork)
    .check_suggested(pkg = "forcats")
    .check_suggested(pkg = "ggridges")
    .check_suggested(pkg = "patchwork")

    if (!"sensitivity_analysis" %in% names(
        MultiAssayExperiment::metadata(expomicset)$differential_analysis
    )) {
        stop("Please run `run_sensitivity_analysis()` first.")
    }

    feature_stability <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis$feature_stability

    if (!stability_metric %in% colnames(feature_stability)) {
        stop(sprintf(
            "Invalid stability_metric: '%s'. Must be one of: %s",
            stability_metric,
            paste(colnames(feature_stability), collapse = ", ")
        ))
    }

    # sensitivity_sum <- feature_stability |>
    #   dplyr::group_by(exp_name) |>
    #   dplyr::summarise(n=dplyr::n()) |>
    #   dplyr::arrange(desc(n))

    sensitivity_sum <- feature_stability |>
        dplyr::group_by(exp_name) |>
        dplyr::summarise(
            n_above = sum(
                .data[[stability_metric]] > stability_score_thresh,
                na.rm = TRUE
            ),
            total = dplyr::n()
        ) |>
        dplyr::arrange(desc(total))

    sensitivity_bar <- sensitivity_sum |>
        ggplot(aes(
            x = n_above,
            y = forcats::fct_reorder(exp_name, n_above),
            fill = exp_name
        )) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_segment(aes(
            x = n_above,
            xend = n_above,
            y = as.numeric(forcats::fct_reorder(exp_name, n_above)) - 0.45,
            yend = as.numeric(forcats::fct_reorder(exp_name, n_above)) + 0.45,
            color = exp_name,
        ), size = 1) +
        scale_fill_tidy_exp() +
        scale_color_tidy_exp() +
        ggpubr::theme_pubr(legend = "none") +
        theme(
            plot.title = element_text(face = "bold.italic"),
            axis.text.y = element_blank(),
            text = element_text(size = 10)
        ) +
        labs(
            title = "",
            y = "",
            x = "No. of Features"
        )

    if (is.null(stability_score_thresh)) {
        stability_score_thresh <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis$score_thresh
    }

    sensitivity_ridgeplot <- feature_stability |>
        dplyr::left_join(sensitivity_sum,
            by = "exp_name"
        ) |>
        mutate(exp_name = paste0(exp_name, ": (", n_above, "/", total, ")")) |>
        ggplot(aes(
            x = .data[[stability_metric]],
            y = forcats::fct_reorder(exp_name, n_above),
            fill = exp_name
        )) +
        ggridges::geom_density_ridges() +
        scale_fill_tidy_exp() +
        theme_minimal() +
        geom_vline(
            xintercept = stability_score_thresh,
            linetype = "dashed",
            color = "grey55"
        ) +
        theme(
            plot.title = element_text(face = "bold.italic"),
            legend.position = "none",
            axis.text.y = element_text(size = 10, color = "black")
        ) +
        labs(
            fill = "Assay",
            x = "Stability Score",
            y = "",
            title = title
        )


    message(
        "Number of Features with ",
        stability_metric,
        " > ",
        stability_score_thresh, ":"
    )
    for (i in seq_len(nrow(sensitivity_sum))) {
        message(
            sensitivity_sum$exp_name[i],
            ": ",
            sensitivity_sum$n_above[i],
            " / ",
            sensitivity_sum$total[i]
        )
    }

    return((sensitivity_ridgeplot | sensitivity_bar) +
        patchwork::plot_layout(
            widths = c(3, 1)
        ))
}
