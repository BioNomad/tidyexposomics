#' Plot Sample or Exposure Clusters
#'
#' Generates a heatmap of either sample clustering or exposure clustering
#' results, split by cluster group.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing clustering
#' results in `metadata(exposomicset)$quality_control`.
#' @param cluster_type A character string specifying which clustering to plot:
#' `"sample"` (from `run_cluster_samples()`) or `"exposure"`
#' (from `run_cluster_exposures()`). Default is `"sample"`.
#' @param exposure_cols A character vector specifying columns from `colData`
#' to include in the heatmap. Default is `NULL`, which includes all `colData`
#' columns for `"sample"` and all clustered exposures for `"exposure"`.
#'
#' @details
#' This function:
#' - Extracts cluster assignments from
#' `metadata(exposomicset)$quality_control`.
#' - Merges cluster labels with `colData(exposomicset)` and the codebook.
#' - Plots a heatmap split by sample groups (columns) or exposure groups
#' (rows), with annotation tiles for exposure category and cluster group.
#'
#' @return A `ComplexHeatmap` plot displaying clustering results.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 30,
#'     return_mae = TRUE
#' )
#'
#' exposure_cols <- c("exposure_pm25", "exposure_no2", "age", "bmi")
#'
#' # determine sample and exposure clusters
#' mae <- mae |>
#'     run_cluster_samples(
#'         exposure_cols = exposure_cols,
#'         clustering_approach = "diana"
#'     ) |>
#'     run_cluster_exposures(
#'         exposure_cols = exposure_cols,
#'         clustering_approach = "diana"
#'     )
#'
#' # plot sample clusters
#' sample_cluster_p <- mae |>
#'     plot_clusters(
#'         cluster_type = "sample",
#'         exposure_cols = exposure_cols
#'     )
#'
#' # plot exposure clusters
#' exposure_cluster_p <- mae |>
#'     plot_clusters(
#'         cluster_type = "exposure",
#'         exposure_cols = exposure_cols
#'     )
#'
#' @importFrom dplyr select mutate across inner_join rename group_by all_of
#' @importFrom dplyr everything
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom purrr pluck
#' @importFrom rlang sym
#' @export
plot_clusters <- function(
  exposomicset,
  cluster_type = c("sample", "exposure"),
  exposure_cols = NULL
) {
    cluster_type <- match.arg(cluster_type)

    # check for suggested packages
    .check_suggested(pkg = "tidyHeatmap")
    .check_suggested(pkg = "circlize")

    cluster_key <- paste0(cluster_type, "_clustering")
    group_col <- paste0(cluster_type, "_group")

    if (!cluster_key %in% names(
        MultiAssayExperiment::metadata(exposomicset)$quality_control
    )) {
        stop("Please run `run_cluster_", cluster_type, "s()` first")
    }

    groups <- MultiAssayExperiment::metadata(exposomicset) |>
        purrr::pluck(
            "quality_control",
            cluster_key,
            paste0(cluster_type, "_groups")
        )

    # sample groups map onto columns (Sample), exposure groups onto rows
    join_col <- if (cluster_type == "sample") "Sample" else "variable"
    group_prefix <- "G"

    if (is.null(exposure_cols)) {
        exposure_cols <- if (cluster_type == "sample") {
            colnames(MultiAssayExperiment::colData(exposomicset))
        } else {
            names(groups)
        }
    }

    exposure_cols <- intersect(
        exposure_cols,
        colnames(MultiAssayExperiment::colData(exposomicset))
    )

    df <- exposomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame() |>
        dplyr::select(dplyr::all_of(exposure_cols)) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
        tibble::rownames_to_column("Sample") |>
        tidyr::pivot_longer(
            -Sample,
            names_to = "variable",
            values_to = "value"
        ) |>
        dplyr::inner_join(
            MultiAssayExperiment::metadata(exposomicset)$codebook |>
                dplyr::select(dplyr::all_of(c("variable", "category"))),
            by = "variable"
        ) |>
        dplyr::rename(Category = category) |>
        dplyr::inner_join(
            data.frame(
                names(groups),
                paste0(group_prefix, as.character(groups))
            ) |>
                stats::setNames(c(join_col, group_col)),
            by = join_col
        )

    df |>
        dplyr::mutate(Category = as.character(Category)) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
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
            !!rlang::sym(group_col),
            palette = tidy_exp_pal[
                seq_len(length(unique(df[[group_col]])))
            ]
        )
}
