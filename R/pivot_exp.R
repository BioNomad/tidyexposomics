#' Pivot a selected omics dataset from a MultiAssayExperiment into tidybulk format
#'
#' Extracts a specified omics dataset from a \code{MultiAssayExperiment},
#' optionally filters by feature (row) names,
#' and returns a tidy tibble. The output includes assay values along with
#' sample metadata and feature metadata.
#'
#' @param exposomicset A \code{MultiAssayExperiment} object containing
#'   one or more omics assays.
#' @param exp_name A character string. The name of the omics dataset
#'   to extract (e.g., "Proteomics").
#' @param features Optional character vector of row (feature) names to retain.
#'   If \code{NULL}, all features are included.
#'
#' @return A tibble in tidy format with one row per feature/sample pair,
#'   including all metadata and a new column \code{exp_name}
#'   indicating the assay source. Assay values are provided in separate columns
#'   named after the assay slot(s).
#'
#' @importFrom dplyr mutate left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom SummarizedExperiment assays rowData colData
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # pivot experiment
#' exp_data <- pivot_exp(
#'     exposomicset = mae,
#'     exp_name = "mRNA",
#'     features = c("feat_42")
#' )
#'
#' @export
pivot_exp <- function(exposomicset,
                      exp_name,
                      features = NULL) {
    # Extract the SummarizedExperiment
    se <- .update_assay_colData(exposomicset, exp_name = exp_name)

    # Optionally subset features
    if (!is.null(features)) {
        se <- se[rownames(se) %in% features, ]
    }

    # Row and column data
    rd <- as.data.frame(rowData(se))
    cd <- as.data.frame(colData(se))

    # Iterate over assays and pivot to long format
    assay_list <- lapply(names(assays(se)), function(aname) {
        a <- assays(se)[[aname]]
        tibble::as_tibble(a, rownames = ".feature") |>
            tidyr::pivot_longer(
                cols = -.feature,
                names_to = ".sample",
                values_to = aname
            )
    })

    # Merge multiple assays by feature/sample pair
    assay_df <- Reduce(
        function(x, y) dplyr::left_join(x, y, by = c(".feature", ".sample")),
        assay_list
    )

    # Add rowData and colData
    assay_df <- assay_df |>
        dplyr::left_join(
            tibble::rownames_to_column(rd, ".feature"),
            by = ".feature"
        ) |>
        dplyr::left_join(
            tibble::rownames_to_column(cd, ".sample"),
            by = ".sample"
        ) |>
        dplyr::mutate(
            exp_name = exp_name,
            .before = .feature
        )

    assay_df
}
