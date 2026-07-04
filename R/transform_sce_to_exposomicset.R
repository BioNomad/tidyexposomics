#' Convert a SingleCellExperiment to an exposomicset
#'
#' Pseudo-bulk aggregates a \code{SingleCellExperiment} to the donor level,
#' splitting by cell type so that each cell type becomes a separate named
#' experiment in the returned \code{MultiAssayExperiment}. Donor-level
#' metadata is extracted directly from \code{colData(sce)}.
#'
#' @param sce A \code{SingleCellExperiment} object with at least one assay
#'   and cell-level metadata columns matching \code{sample_col} and
#'   \code{cluster_col}.
#' @param codebook A data frame with columns \code{variable} and
#'   \code{category} describing the donor-level metadata columns selected
#'   via \code{sample_meta_cols}.
#' @param sample_col Character. Column in \code{colData(sce)} identifying
#'   the donor or sample each cell belongs to. Default \code{"sample"}.
#' @param cluster_col Character. Column in \code{colData(sce)} identifying
#'   cell type or cluster labels. Each unique value becomes a separate
#'   experiment in the MAE. Default \code{"label"}.
#' @param sample_meta_cols Character vector of column names in
#'   \code{colData(sce)} to retain as donor-level metadata. These should be
#'   donor-level variables (e.g. age, sex, smoking status) rather than
#'   cell-level variables (e.g. n_genes_detected). The function deduplicates
#'   to one row per donor after selecting these columns. If \code{NULL},
#'   only \code{sample_col} is used and \code{colData} will be minimal.
#' @param assay_name Character. Assay to aggregate. Default \code{"counts"}.
#' @param min_cells_per_donor Integer. Minimum number of cells a donor must
#'   contribute to a given cell type to be included in that cell type's
#'   pseudo-bulk matrix. Donors below this are dropped from that experiment
#'   only. Default \code{10}.
#' @param min_donor_frac Numeric between 0 and 1. Minimum fraction of donors
#'   that must pass \code{min_cells_per_donor} for a cell type to be retained
#'   as an experiment. Cell types below this threshold are dropped with a
#'   warning. Default \code{0.8}.
#'
#' @return A \code{MultiAssayExperiment} exposomicset with one pseudo-bulk
#'   experiment per retained cell type, ready for
#'   \code{run_exposure_omics_association()}, \code{run_correlation()},
#'   \code{run_association()}, and all other tidyexposomics functions.
#'
#' @details
#' Pseudo-bulk aggregation sums raw counts across all cells belonging to the
#' same donor and cell type.
#'
#'
#'
#' @examples
#' \dontrun{
#' mae <- transform_sce_to_exposomicset(
#'     sce              = sce,
#'     codebook         = my_codebook,
#'     sample_col       = "donor_id",
#'     cluster_col      = "cell_type",
#'     sample_meta_cols = c("age", "sex", "smoking_status", "BMI")
#' )
#'
#' mae <- mae |>
#'     run_exposure_omics_association(
#'         exposures  = c("pm25", "no2"),
#'         covariates = c("age", "sex")
#'     )
#' }
#'
#' @export
transform_sce_to_exposomicset <- function(
  sce,
  codebook,
  sample_col = "sample",
  cluster_col = "label",
  sample_meta_cols = NULL,
  assay_name = "counts",
  min_cells_per_donor = 10,
  min_donor_frac = 0.8
) {
    # Validate
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment.")
    }
    if (!sample_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop(sprintf("Column '%s' not found in colData(sce).", sample_col))
    }
    if (!cluster_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop(sprintf("Column '%s' not found in colData(sce).", cluster_col))
    }
    if (!is.null(sample_meta_cols)) {
        missing_cols <- setdiff(
            sample_meta_cols,
            colnames(SummarizedExperiment::colData(sce))
        )
        if (length(missing_cols) > 0) {
            stop(
                "These sample_meta_cols were not found in colData(sce): ",
                paste(missing_cols, collapse = ", ")
            )
        }
    }

    # Extract donor-level metadata from colData
    cols_to_keep <- c(sample_col, sample_meta_cols)

    sample_meta <- SummarizedExperiment::colData(sce) |>
        as.data.frame() |>
        dplyr::select(dplyr::all_of(cols_to_keep)) |>
        dplyr::distinct(dplyr::pick(dplyr::all_of(sample_col)), .keep_all = TRUE) |>
        tibble::remove_rownames() |>
        tibble::column_to_rownames(sample_col)

    n_donors <- nrow(sample_meta)
    cell_types <- unique(SummarizedExperiment::colData(sce)[[cluster_col]])

    message(sprintf(
        "Found %d donors and %d cell types. Aggregating per cell type.",
        n_donors, length(cell_types)
    ))


    # Donor x cell type cell counts for filtering
    cell_meta <- data.frame(
        sample = SummarizedExperiment::colData(sce)[[sample_col]],
        celltype = SummarizedExperiment::colData(sce)[[cluster_col]],
        stringsAsFactors = FALSE
    )

    donor_counts <- cell_meta |>
        dplyr::count(sample, celltype)

    # pull the full assay matrix once — features x cells
    assay_mat <- SummarizedExperiment::assay(sce, assay_name)

    # Aggregate per cell type
    omics_list <- list()
    row_data_list <- list()

    for (ct in cell_types) {
        eligible <- donor_counts |>
            dplyr::filter(celltype == ct, n >= min_cells_per_donor) |>
            dplyr::pull(sample)

        frac_eligible <- length(eligible) / n_donors

        if (frac_eligible < min_donor_frac) {
            warning(sprintf(
                "Dropping '%s': only %.0f%% of donors have >= %d cells (need %.0f%%).",
                ct, frac_eligible * 100, min_cells_per_donor, min_donor_frac * 100
            ))
            next
        }

        # cells belonging to this cell type from eligible donors
        cells_keep <- cell_meta$celltype == ct & cell_meta$sample %in% eligible
        mat_ct <- assay_mat[, cells_keep, drop = FALSE]
        donors_ct <- cell_meta$sample[cells_keep]

        # sum counts per donor, vapply over unique donors, rowSums per group
        pb_mat <- vapply(eligible, function(d) {
          cols <- donors_ct == d
          if (sum(cols) == 1L) {
            as.numeric(mat_ct[, cols, drop = TRUE])
          } else {
            as.numeric(Matrix::rowSums(mat_ct[, cols, drop = FALSE]))
          }
        }, FUN.VALUE = numeric(nrow(assay_mat)))

        # ensure matrix with correct dimnames
        pb_mat <- as.matrix(pb_mat)
        rownames(pb_mat) <- rownames(assay_mat)
        colnames(pb_mat) <- eligible

        # align columns to sample_meta rownames
        pb_mat <- pb_mat[, intersect(colnames(pb_mat), rownames(sample_meta)), drop = FALSE]

        ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
        omics_list[[ct_safe]] <- pb_mat

        rd <- SummarizedExperiment::rowData(sce) |> as.data.frame()
        row_data_list[[ct_safe]] <- S4Vectors::DataFrame(
            rd[rownames(pb_mat), , drop = FALSE]
        )

        message(sprintf(
            "  '%s': %d genes x %d donors.",
            ct_safe, nrow(pb_mat), ncol(pb_mat)
        ))
    }

    if (length(omics_list) == 0) {
        stop(
            "No cell types passed the filters. ",
            "Try lowering min_donor_frac or min_cells_per_donor."
        )
    }

    message(sprintf(
        "%d / %d cell types retained.",
        length(omics_list), length(cell_types)
    ))

    # Build MAE
    create_exposomicset(
        codebook = codebook,
        exposure = sample_meta,
        omics    = omics_list,
        row_data = row_data_list
    )
}
