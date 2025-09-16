#' Plot Sample Clusters
#'
#' Generates a heatmap of sample clustering results and
#' summarizes sample group assignments.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing sample
#'  clustering results
#' in `metadata(exposomicset)$sample_clustering`.
#' @param exposure_cols A character vector specifying columns from `colData`
#' to include
#' in the summary. Default is `NULL`, which includes all available columns.
#'
#' @details
#' This function:
#' - Extracts sample cluster assignments from
#' `metadata(exposomicset)$sample_clustering`.
#' - Merges cluster labels with `colData(exposomicset)`.
#' - Plots the heatmap stored in
#' `metadata(exposomicset)$sample_clustering$heatmap`.
#'
#' @return A `ComplexHeatmap` plot displaying sample clustering results.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 30,
#'     return_mae = TRUE
#' )
#'
#' # determine sample clusters
#' mae <- run_cluster_samples(
#'     exposomicset = mae,
#'     exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
#'     clustering_approach = "diana"
#' )
#'
#' # plot sample clusters
#' sample_cluster_p <- mae |>
#'     plot_sample_clusters(
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' @importFrom dplyr select mutate_all inner_join rename group_by mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom purrr pluck
#' @export
plot_sample_clusters <- function(
    exposomicset,
    exposure_cols = NULL) {
    # check for suggested packages
    .check_suggested(pkg = "tidyHeatmap")
    .check_suggested(pkg = "circlize")
    if (!"sample_clustering" %in% names(
        MultiAssayExperiment::metadata(exposomicset)$quality_control
    )) {
        stop("Please run `run_cluster_samples()` first")
    }

    if (is.null(exposure_cols)) {
        exposure_cols <- colnames(exposomicset)
    } else {
        exposure_cols <- exposure_cols[exposure_cols %in% colnames(
            MultiAssayExperiment::colData(exposomicset)
        )]
    }

    df <- exposomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame() |>
        dplyr::select(exposure_cols) |>
        dplyr::mutate_all(~ as.numeric(.)) |>
        tibble::rownames_to_column("Sample") |>
        dplyr::inner_join(
            data.frame(
                Sample = names(
                    MultiAssayExperiment::metadata(exposomicset) |>
                        purrr::pluck(
                            "quality_control",
                            "sample_clustering",
                            "sample_groups"
                        )
                ),
                sample_group = paste0(
                    "Group_",
                    as.character(
                        MultiAssayExperiment::metadata(exposomicset) |>
                            purrr::pluck(
                                "quality_control",
                                "sample_clustering",
                                "sample_groups"
                            )
                    )
                )
            ),
            by = c("Sample" = "Sample")
        ) |>
        tidyr::pivot_longer(-c(Sample, sample_group),
            names_to = "variable",
            values_to = "value"
        ) |>
        dplyr::inner_join(
            MultiAssayExperiment::metadata(exposomicset)$codebook,
            by = "variable"
        ) |>
        dplyr::rename(Category = category)

    df |>
        mutate(Category = as.character(Category)) |>
        dplyr::group_by(sample_group) |>
        tidyHeatmap::heatmap(
            variable,
            Sample,
            value,
            scale = "row",
            palette_value = circlize::colorRamp2(
                c(-2, 0, 2),
                c("#006666", "white", "#8E0152")
            )
        ) |>
        tidyHeatmap::annotation_tile(
            Category,
            palette = rev(tidy_exp_pal)[
                seq_len(length(unique(df$Category)))
            ]
        ) |>
        tidyHeatmap::annotation_tile(
            sample_group,
            palette = tidy_exp_pal[
                seq_len(length(unique(df$sample_group)))
            ]
        )
}
