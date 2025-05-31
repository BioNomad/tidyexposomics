#' Filter Sample Outliers
#'
#' Removes sample outliers from a `MultiAssayExperiment` object based on PCA analysis.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param outliers An optional character vector specifying sample names to be removed.
#' If `NULL`, the function uses outliers identified in `metadata(expomicset)$pca$outliers`. Default is `NULL`.
#'
#' @details
#' The function checks for the presence of PCA results in `metadata(expomicset)`. If `outliers` is not provided,
#' it retrieves precomputed outliers from `metadata(expomicset)$pca$outliers`. The identified samples are removed
#' from the dataset.
#'
#' @return A `MultiAssayExperiment` object with the specified outliers removed.
#'
#' @examples
#' \dontrun{
#' filtered_expomicset <- filter_sample_outliers(
#'   expomicset = my_expomicset
#' )
#' }
#'
#' @export
filter_sample_outliers <- function(
    expomicset,
    outliers=NULL
){

  # Check if the input is a MultiAssayExperiment object
  if(!"pca" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("PCA not performed on the data. Please run 'pca_analysis' first.")
  }

  # Check if outliers are provided
  if(is.null(outliers)){
    outliers <-  MultiAssayExperiment::metadata(expomicset)$pca$outliers
  } else{
    outliers <- outliers
  }

  message("Removing outliers: ", paste(outliers, collapse=", "))
  # Subset by ColData to properly filter outliers
  expomicset <- MultiAssayExperiment::subsetByColData(expomicset,
                                                      !(rownames(SummarizedExperiment::colData(expomicset)) %in% outliers))

  # Add analysis steps taken to metadata
  MultiAssayExperiment::metadata(expomicset)$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$steps,
    "filter_sample_outliers"
  )

  # expomicset <- expomicset[,!colnames(expomicset) %in% outliers]
  return(expomicset)
}
