#' Pivot a selected omics dataset from a MultiAssayExperiment into tidybulk format
#'
#' Extracts a specified omics dataset from a \code{MultiAssayExperiment},
#' optionally filters by feature (row) names,
#' and returns a tidy tibble. The output includes assay values along with
#' sample metadata and feature metadata.
#'
#' @param exposomicset A \code{MultiAssayExperiment} object containing
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
#' @importFrom dplyr mutate
#' @importMethodsFrom tidybulk tidybulk
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
#' @export
pivot_exp <- function(
    exposomicset,
    omics_name,
    features = NULL) {
    # Extract the SummarizedExperiment
    se <- .update_assay_colData(
        exposomicset,
        exp_name = omics_name
    )

    # Optionally subset features
    if (!is.null(features)) {
        se <- se[rownames(se) %in% features, ]
    }

    # Apply tidybulk (uses assay, colData, rowData internally)
    tidybulk(se) |>
        mutate(exp_name = omics_name, .before = .feature)
}
