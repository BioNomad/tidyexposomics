#' Plot Missing Data Across Exposure and Omic Layers
#'
#' Visualizes missing data patterns in a `MultiAssayExperiment` object using
#' summary bar plots or feature-level lollipop plots.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure
#' and omics assays. Missing data is inferred directly from the assays.
#' @param threshold Numeric. The percentage threshold (0–100) above which
#' features are counted as missing in the summary plot. Default is `5`.
#' @param plot_type Character. Type of plot to generate. Either `"summary"`
#' for a bar plot showing number of features above the missing threshold,
#'  or `"lollipop"` for a per-feature lollipop plot with layer annotations.
#'   Default is `"summary"`.
#' @param layers Optional character vector. If specified, filters the plot
#' to include only selected layers (e.g., `"Exposure"`, `"Transcriptome"`).
#'
#' @details
#' The function calculates missing data per feature (or variable) across all
#'  assays (including exposure variables) and generates:
#'
#' - **Summary plot (`plot_type = "summary`)**: A bar plot showing the number
#' of variables in each assay exceeding the specified missingness threshold.
#' - **Lollipop plot (`plot_type = "lollipop`)**: A feature-level plot where
#' each feature's percent missingness is shown, along with a color-coded tile
#' on the side indicating its layer of origin.
#'
#' The tile colors in the lollipop plot match the experiment colors used in
#' other visualizations (e.g., via `scale_color_tidy_exp()`).
#'
#' @return A `ggplot` or `patchwork` object depending on the selected
#' `plot_type`.
#'
#' @import ggplot2
#'
#' @examples
#'
#' #' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # Introduce some missingness
#' MultiAssayExperiment::colData(mae)$exposure_pm25[sample(1:20, 5)] <- NA
#'
#' # Summary bar plot of missing data
#' summary_p <- plot_missing(
#'     mae,
#'     threshold = 10,
#'     plot_type = "summary"
#' )
#'
#' # Lollipop plot for all features with any missingness
#' lollipop_p <- plot_missing(
#'     mae,
#'     plot_type = "lollipop"
#' )
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_point geom_tile
#'   labs theme_bw theme element_text theme_void
#' @importFrom MultiAssayExperiment colData experiments
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr bind_rows filter group_by reframe mutate arrange
#'   reframe select ungroup
#' @importFrom purrr map
#' @importFrom naniar miss_var_summary
#' @importFrom ggpubr theme_pubr
#'
#' @export
plot_missing <- function(
    expomicset,
    threshold = 5,
    plot_type = c("summary", "lollipop"),
    layers = NULL) {
    # require(ggplot2)
    # require(patchwork)
    .check_suggested(pkg = "forcats")
    .check_suggested(pkg = "patchwork")

    plot_type <- match.arg(plot_type)
    # grab exposure data
    exposure <- MultiAssayExperiment::colData(expomicset) |>
        as.data.frame()

    # grab omics assay data
    omics <- lapply(
        names(MultiAssayExperiment::experiments(expomicset)),
        function(exp_name) {
            exp <- expomicset[[exp_name]]
            assay <- SummarizedExperiment::assay(exp) |>
                t() |>
                as.data.frame()
        }
    ) |>
        setNames(names(MultiAssayExperiment::experiments(expomicset)))

    # get the number of NA values per exposure variable and omic feature
    na_df <- c(
        list(Exposure = exposure),
        omics
    ) |>
        purrr::map(~ .x |>
            naniar::miss_var_summary()) |>
        dplyr::bind_rows(.id = "exp_name")

    if (!is.null(layers)) {
        na_df <- na_df |>
            dplyr::filter(exp_name %in% layers)
    }

    if (plot_type == "summary") {
        plot_data <- na_df |>
            dplyr::group_by(exp_name) |>
            dplyr::reframe(
                missingness = sum(pct_miss > threshold),
                .groups = "drop"
            ) |>
            dplyr::mutate(exp_label = paste0(exp_name, " (", missingness, ")"))

        p <- plot_data |>
            ggplot(aes(
                y = forcats::fct_reorder(exp_label, missingness),
                x = missingness,
                fill = exp_name
            )) +
            geom_bar(stat = "identity", alpha = 0.7) +
            geom_segment(aes(
                x = missingness,
                xend = missingness,
                y = as.numeric(forcats::fct_reorder(exp_label, missingness)) - 0.45,
                yend = as.numeric(forcats::fct_reorder(exp_label, missingness)) + 0.45,
                color = exp_name
            ), size = 1) +
            scale_fill_tidy_exp() +
            scale_color_tidy_exp() +
            labs(
                title = paste("No. of Features Over ", threshold, "% Threshold"),
                y = NULL,
                x = "No. of Features"
            ) +
            ggpubr::theme_pubr(legend = "none")
    } else if (plot_type == "lollipop") {
        na_df_filtered <- na_df |>
            dplyr::filter(pct_miss > 0) |>
            dplyr::arrange(pct_miss) |>
            dplyr::mutate(var_label = factor(variable, levels = unique(variable)))

        # main lollipop plot
        p_lollipop <- na_df_filtered |>
            ggplot(aes(
                x = pct_miss,
                y = var_label
            )) +
            geom_segment(
                aes(
                    x = 0,
                    xend = pct_miss,
                    yend = var_label
                ),
                color = "midnightblue"
            ) +
            geom_point(color = "midnightblue") +
            theme_bw(base_size = 10) +
            labs(
                title = "Missingness per Feature",
                x = "% Missing",
                y = NULL
            ) +
            theme(
                axis.text.y = element_text(),
                plot.title = element_text(face = "bold.italic"),
                legend.position = "right"
            )

        # side tile plot (just layer colors)
        p_tile <- na_df_filtered |>
            ggplot(aes(
                x = 1,
                y = var_label,
                fill = exp_name
            )) +
            geom_tile() +
            scale_fill_tidy_exp() +
            theme_void() +
            theme(legend.position = "right") +
            labs(fill = "Experiment")

        # combine plots with patchwork
        p <- p_lollipop + p_tile + patchwork::plot_layout(
            widths = c(1, 0.05),
            guides = "collect"
        )
    }

    p
}
