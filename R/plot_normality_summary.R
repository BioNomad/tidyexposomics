#' Plot Normality Summary of Exposure Variables
#'
#' Generates a bar plot summarizing the number of exposure variables that
#' pass or fail normality
#' tests (e.g., Shapiro-Wilk) before or after transformation.
#'
#' @param expomicset A `MultiAssayExperiment` object with quality control
#'  metadata.
#' @param transformed Logical; if `TRUE`, use results after transformation.
#' Default is `FALSE`.
#'
#' @return A `ggplot` object summarizing the number of exposures
#' classified as normal or not normal.
#'
#' @details
#' This function assumes that `run_normality_check()` has been executed and
#'  that the results are
#' stored in `metadata(expomicset)$quality_control$normality`.
#'  If `transformed = TRUE`, the function will
#' instead plot the transformation summary stored in `metadata(expomicset)$quality_control$transformation$norm_summary`,
#' which is populated by `transform_exposure()`.
#'
#' The plot includes both bar heights and overlaid line segments to
#' reinforce the counts.
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # Test for normality
#' mae <- mae |>
#'     run_normality_check() |>
#'     transform_exposure(exposure_cols = c("age", "bmi", "exposure_pm25"))
#'
#' # plot the normality summary
#' norm_p <- mae |>
#'     plot_normality_summary()
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment labs scale_fill_manual scale_color_manual theme
#' @importFrom ggpubr theme_pubr get_palette
#' @importFrom purrr pluck
#' @export
plot_normality_summary <- function(
    expomicset,
    transformed = FALSE) {
    # require(ggplot2)

    # Check if "normality" is a name in metadata
    if (!("normality" %in% names(
        MultiAssayExperiment::metadata(expomicset)$quality_control
    ))) {
        stop("Please run `run_normality_check() first.`")
    }

    # Check if "transformation" is a name in metadata
    if (transformed) {
        if (!("transformation" %in% names(
            MultiAssayExperiment::metadata(expomicset)$quality_control
        ))) {
            stop("Please run `transform_exposure() first.`")
        }
    }

    # Plot normality results
    if (transformed) {
        norm_plot <- MultiAssayExperiment::metadata(expomicset)$quality_control$transformation$norm_summary |>
            ggplot(aes(
                x = var,
                y = value,
                fill = var
            )) +
            geom_bar(stat = "identity", alpha = 0.5) +
            geom_segment(
                aes(
                    x = as.numeric(as.factor(var)) - 0.45,
                    xend = as.numeric(as.factor(var)) + 0.45,
                    y = value,
                    yend = value,
                    color = var
                ),
                linewidth = 1
            ) +
            ggpubr::theme_pubr(legend = "right") +
            scale_fill_manual(
                values = ggpubr::get_palette("uchicago", k = 2)[c(2, 1)]
            ) +
            scale_color_manual(
                values = ggpubr::get_palette("uchicago", k = 2)[c(2, 1)]
            ) +
            guides(color = FALSE) +
            theme(
                plot.title = element_text(face = "bold.italic"),
                plot.subtitle = element_text(face = "italic")
            ) +
            labs(
                x = "",
                y = "No. of Exposures",
                fill = "",
                title = "Normality of Exposure Variables",
                subtitle = "Shapiro-Wilk Test"
            )
    } else {
        norm_plot <- MultiAssayExperiment::metadata(expomicset) |>
            purrr::pluck(
                "quality_control",
                "normality",
                "norm_summary"
            ) |>
            ggplot(aes(
                x = var,
                y = value,
                fill = var
            )) +
            geom_bar(stat = "identity", alpha = 0.5) +
            geom_segment(
                aes(
                    x = as.numeric(as.factor(var)) - 0.45,
                    xend = as.numeric(as.factor(var)) + 0.45,
                    y = value,
                    yend = value,
                    color = var
                ),
                linewidth = 1
            ) +
            ggpubr::theme_pubr(legend = "right") +
            scale_fill_tidy_exp() +
            scale_color_tidy_exp(guide = "none") +
            # ggsci::scale_fill_lancet() +
            # ggsci::scale_color_lancet(guide = FALSE) +
            theme(
                plot.title = element_text(face = "bold.italic"),
                plot.subtitle = element_text(face = "italic")
            ) +
            labs(
                x = "",
                y = "No. of Exposures",
                fill = "",
                title = "Normality of Exposure Variables",
                subtitle = "Shapiro-Wilk Test"
            )
    }

    return(norm_plot)
}
