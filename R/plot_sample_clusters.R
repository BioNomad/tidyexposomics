#' Plot Sample Clusters
#'
#' Generates a heatmap of sample clustering results and summarizes sample group assignments.
#'
#' @param expomicset A `MultiAssayExperiment` object containing sample clustering results 
#' in `metadata(expomicset)$sample_clustering`.
#' @param cols_of_interest A character vector specifying columns from `colData` to include 
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
    cols_of_interest=NULL
){
  
  if(!"sample_clustering" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `cluster_samples()` first")
  }
  
  if(is.null(cols_of_interest)){
    cols_of_interest <- colnames(expomicset)
  } else{
    cols_of_interest <- cols_of_interest[cols_of_interest %in% colnames(MultiAssayExperiment::colData(expomicset))]
  }
  
  # add in sample group to the colData
  meta <- expomicset |> 
    MultiAssayExperiment::colData() |> 
    as.data.frame() |> 
    tibble::rownames_to_column("id_to_map") |>
    dplyr::left_join(
      data.frame(
        id_to_map = names(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups),
        cluster = as.vector(MultiAssayExperiment::metadata(expomicset)$sample_clustering$sample_groups)
      ),
      by="id_to_map"
    ) |> 
    dplyr::select(dplyr::all_of(c(cols_of_interest,"cluster"))) |> 
    dplyr::mutate_at(dplyr::vars(cols_of_interest), as.numeric)
  
  ComplexHeatmap::draw(MultiAssayExperiment::metadata(expomicset)$sample_clustering$heatmap)
  
}





