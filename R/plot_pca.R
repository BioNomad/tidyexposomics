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
#' @import ggplot2
#' @importFrom ggpubr theme_pubr
#' @importFrom factoextra fviz_eig
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot_pca <- function(
    expomicset,
    feature_col = "#00a9b2",
    sample_col = "#8a4f77",
    sample_outlier_col = "firebrick") {
    .check_suggested(pkg = "patchwork")

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
    pca_plot_feature <- pca_feature |>
        purrr::pluck("x") |>
        as.data.frame() |>
        mutate(category = dat$category) |>
        ggplot(aes(
            x = PC1,
            y = PC2,
            color = category
        )) +
        geom_point() +
        ggpubr::theme_pubr(legend = "right") +
        scale_color_tidy_exp() +
        labs(
            title = "PCA of Feature Space",
            color = "",
            x = sprintf(
                "PC1 (%s%%)",
                round(
                    (summary(pca_feature)$importance[2, seq_len(2)] * 100)[1],
                    digits = 1
                )
            ),
            y = sprintf(
                "PC2 (%s%%)",
                round(
                    (summary(pca_feature)$importance[2, seq_len(2)] * 100)[2],
                    digits = 1
                )
            )
        ) +
        theme(
            plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic")
        )
    # pca_plot_feature <- autoplot(pca_feature,
    #     data = dat,
    #     colour = "category"
    # ) +
    #     ggpubr::theme_pubr(legend = "right") +
    #     ggsci::scale_color_aaas() +
    #     labs(title = "PCA of Feature Space", color = "") +
    #     theme(
    #         plot.title = element_text(face = "bold.italic"),
    #         plot.subtitle = element_text(face = "italic")
    #     )

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
    pca_plot_sample <- pca_sample |>
        purrr::pluck("x") |>
        as.data.frame() |>
        (\(df){
            df$id <- rownames(df)
            df
        })() |>
        mutate(label = ifelse(id %in% outlier_samples, id, NA)) |>
        ggplot(aes(
            x = PC1,
            y = PC2
        )) +
        geom_point(color = sample_col) +
        ggrepel::geom_text_repel(
            aes(
                x = PC1,
                y = PC2,
                label = label
            ),
            color = sample_outlier_col
        ) +
        ggpubr::theme_pubr() +
        labs(
            title = "PCA of Sample Space",
            x = sprintf(
                "PC1 (%s%%)",
                round(
                    (summary(pca_sample)$importance[2, seq_len(2)] * 100)[1],
                    digits = 1
                )
            ),
            y = sprintf(
                "PC2 (%s%%)",
                round(
                    (summary(pca_sample)$importance[2, seq_len(2)] * 100)[2],
                    digits = 1
                )
            )
        ) +
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
