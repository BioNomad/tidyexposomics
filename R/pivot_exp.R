#' Pivot a selected omics dataset from a MultiAssayExperiment into tidybulk format
#'
#' Extracts a specified omics dataset from a \code{MultiAssayExperiment},
#' optionally filters by feature (row) names,
#' and returns a tidy tibble in the structure of
#' \code{tidybulk::tidybulk()}. The output includes assay values along with
#' sample metadata
#' (from \code{colData}) and feature metadata (from \code{rowData}).
#'
#' @param expomicset A \code{MultiAssayExperiment} object containing
#' omics assays.
#' @param omics_name A character string. The name of the omics dataset
#' to extract (e.g., "Proteomics").
#' @param features Optional character vector of row (feature) names to retain.
#'  If \code{NULL}, all features are included.
#'
#' @return A tibble in tidybulk format with one row per feature/sample pair,
#'  including all metadata and a new column \code{exp_name}
#'  indicating the assay source.
#'
#' @importFrom tidybulk tidybulk
#' @importFrom dplyr mutate
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # pivot experiment
#' exp_data <- mae |>
#'     pivot_exp(
#'         omics_name = "mRNA",
#'         features = "feat_42"
#'     )
#'
#' @export
pivot_exp <- function(
    expomicset,
    omics_name,
    features = NULL) {
    # Extract the SummarizedExperiment
    se <- .update_assay_colData(
        expomicset,
        exp_name = omics_name
    )

    # Optionally subset features
    if (!is.null(features)) {
        se <- se[rownames(se) %in% features, ]
    }

    # Apply tidybulk (uses assay, colData, rowData internally)
    tidybulk::tidybulk(se) |>
        mutate(exp_name = omics_name, .before = .feature)
}
