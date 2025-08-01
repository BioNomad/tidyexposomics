#' Plot PCA Results for Features and Samples
#'
#' Generates PCA plots for both feature space and sample space,
#' including scatter plots and scree plots.
#'
#' @param expomicset A `MultiAssayExperiment` object containing PCA results
#' in `metadata(expomicset)$pca`.
#' @param feature_col A character string specifying the color for the
#'  feature scree plot.
#' Default is `"#00a9b2"`.
#' @param sample_col A character string specifying the color for the
#'  sample scree plot.
#' Default is `"#8a4f77"`.
#' @param sample_outlier_col A character string specifying the color
#' for sample outlier labels.
#' Default is `"firebrick"`.
#'
#' @details
#' This function creates four PCA visualizations:
#' - **Feature Space PCA Plot**: Colored by category (e.g., omics, exposure).
#' - **Feature Scree Plot**: Displays the variance explained by each
#' principal component.
#' - **Sample Space PCA Plot**: Highlights outlier samples.
#' - **Sample Scree Plot**: Displays variance explained in the sample PCA.
#'
#' Outliers are labeled based on `metadata(expomicset)$pca$outliers`.
#'
#' @return A combined `ggplot` object containing the four PCA plots.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run pca
#' mae <- mae |>
#'     run_pca()
#'
#' # create the pca plot
#' pca_p <- mae |>
#'     plot_pca()
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom dplyr mutate case_when
#' @importFrom ggplot2 autoplot
#' @import ggfortify
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsci scale_color_aaas
#' @importFrom factoextra fviz_eig
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_pca <- function(
    expomicset,
    feature_col = "#00a9b2",
    sample_col = "#8a4f77",
    sample_outlier_col = "firebrick") {
    # require(ggplot2)
    # require(ggfortify)

    # Check if the required metadata is present
    if (is.null(MultiAssayExperiment::metadata(expomicset)$quality_control$pca)) {
        stop("Please run `run_pca` first.")
    }

    # grab data
    dat <- MultiAssayExperiment::metadata(expomicset)$quality_control$pca$pca_df |>
        as.data.frame()
    pca_feature <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "quality_control",
            "pca",
            "pca_feature"
        )
    pca_sample <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "quality_control",
            "pca",
            "pca_sample"
        )

    # capitalize exposure for plot
    dat <- dat |>
        dplyr::mutate(category = dplyr::case_when(
            category == "exposure" ~ "Exposure",
            .default = category
        ))


    # PCA Feature Scatter Plot
    pca_plot_feature <- autoplot(pca_feature,
        data = dat,
        colour = "category"
    ) +
        ggpubr::theme_pubr(legend = "right") +
        ggsci::scale_color_aaas() +
        labs(title = "PCA of Feature Space", color = "") +
        theme(
            plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic")
        )

    # PCA Feature Scree Plot
    scree_feature <- factoextra::fviz_eig(
        pca_feature,
        barfill = feature_col,
        barcolor = feature_col,
        main = "Scree Plot of Feature Space"
    ) +
        theme(
            plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic")
        )

    # Define outliers
    outlier_samples <- MultiAssayExperiment::metadata(expomicset) |>
        purrr::pluck(
            "quality_control",
            "pca",
            "outliers"
        )

    # PCA sample scatter plot
    pca_plot_sample <- autoplot(pca_sample, colour = sample_col) +
        ggrepel::geom_text_repel(
            aes(
                x = PC1,
                y = PC2,
                label = ifelse(rownames %in% outlier_samples, rownames, NA)
            ),
            color = sample_outlier_col
        ) +
        ggpubr::theme_pubr() +
        labs(title = "PCA of Sample Space") +
        theme(
            plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic")
        )


    # Sample Scree Plot
    scree_sample <- factoextra::fviz_eig(
        pca_sample,
        barfill = sample_col,
        barcolor = sample_col,
        main = "Scree Plot of Sample Space"
    ) +
        theme(
            plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic")
        )

    # Combine plots
    combined_plot <- patchwork::wrap_plots(
        pca_plot_feature, scree_feature,
        pca_plot_sample, scree_sample
    )

    return(combined_plot)
}
