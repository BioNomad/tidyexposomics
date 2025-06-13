#' Plot Sample Clusters
#'
#' Generates a heatmap of sample clustering results and summarizes sample group assignments.
#'
#' @param expomicset A `MultiAssayExperiment` object containing sample clustering results
#' in `metadata(expomicset)$sample_clustering`.
#' @param exposure_cols A character vector specifying columns from `colData` to include
#' in the summary. Default is `NULL`, which includes all available columns.
#'
#' @details
#' This function:
#' - Extracts sample cluster assignments from `metadata(expomicset)$sample_clustering`.
#' - Merges cluster labels with `colData(expomicset)`.
#' - Plots the heatmap stored in `metadata(expomicset)$sample_clustering$heatmap`.
#'
#' @return A `ComplexHeatmap` plot displaying sample clustering results.
#'
#' @examples
#' \dontrun{
#' plot_sample_clusters(expom)
#' }
#'
#' @export
plot_sample_clusters <- function(
    expomicset,
    exposure_cols=NULL
){

  if(!"sample_clustering" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `cluster_samples()` first")
  }

  if(is.null(exposure_cols)){
    exposure_cols <- colnames(expomicset)
  } else{
    exposure_cols <- exposure_cols[exposure_cols %in% colnames(MultiAssayExperiment::colData(expomicset))]
  }

  # add in sample group to the colData
  # meta <- expomicset |>
  #   MultiAssayExperiment::colData() |>
  #   as.data.frame() |>
  #   tibble::rownames_to_column("id_to_map") |>
  #   dplyr::left_join(
  #     data.frame(
  #       id_to_map = names(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups),
  #       cluster = as.vector(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups)
  #     ),
  #     by="id_to_map"
  #   ) |>
  #   dplyr::select(dplyr::all_of(c(exposure_cols,"cluster"))) |>
  #   dplyr::mutate_at(dplyr::vars(exposure_cols), as.numeric)
  #
  # ComplexHeatmap::draw(MultiAssayExperiment::metadata(expomicset)$sample_clustering$heatmap)

  # df <-  a |>
  #   pivot_sample() |>
  #   as.data.frame() |>
  #   dplyr::select(all_of(c(exposure_cols,".sample",sample_group)) |>
  #   mutate(across(exposure_cols,~as.numeric(.))) |>
  #     pivot_longer(-c(.sample,sample_group),
  #                  names_to = "variable",
  #                  values_to = "value") |>
  #     inner_join(
  #       MultiAssayExperiment::metadata(expomicset)$var_info,
  #       by="variable")


  df <- expomicset |>
    MultiAssayExperiment::colData() |>
    as.data.frame() |>
    dplyr::select(exposure_cols) |>
    dplyr::mutate_all(~as.numeric(.)) |>
    tibble::rownames_to_column("Sample") |>
    dplyr::inner_join(
      data.frame(
        Sample = names(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups),
        sample_group = paste0("Group_",
                             as.character(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups))
      ),
      by=c("Sample"="Sample")
    ) |>
    tidyr::pivot_longer(-c(Sample,sample_group),
                 names_to = "variable",
                 values_to = "value") |>
    dplyr::inner_join(
      MultiAssayExperiment::metadata(expomicset)$var_info,
      by="variable") |>
    dplyr::rename(Category=category)


  df |>
    mutate(Category=as.character(Category)) |>
    dplyr::group_by(sample_group) |>
    tidyHeatmap::heatmap(variable,
                         Sample,
                         value,
                         scale = "row",
                         palette_value = circlize::colorRamp2(
                           c(-2, 0, 2),
                           c("#006666", "white", "#8E0152"))) |>
    tidyHeatmap::annotation_tile(
      Category,
      palette = rev(tidy_exp_pal)[1:length(unique(df$Category))]) |>
    tidyHeatmap::annotation_tile(
      sample_group,
      palette = tidy_exp_pal[1:length(unique(df$sample_group))])
}





