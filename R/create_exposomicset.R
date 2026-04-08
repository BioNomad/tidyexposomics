#' Create an Exposomicset Object
#'
#' Constructs a `MultiAssayExperiment` object from exposure data and optionally
#' omics datasets, ensuring proper formatting and alignment of samples and
#' features. For epidemiology-only workflows, omics data can be omitted.
#'
#' @param codebook A data frame containing variable information metadata.
#' @param exposure A data frame containing exposure data,
#' with rows as samples and columns as variables.
#' @param omics An optional list of matrices or a single matrix representing
#' omics data. Each matrix should have samples as columns and features as rows.
#' If `NULL`, creates an exposure-only exposomicset. Default is `NULL`.
#' @param row_data An optional list of `DataFrame` objects providing
#' feature metadata for each omics dataset. If `NULL`,
#' row metadata is generated automatically. Default is `NULL`.
#'
#' @details
#' The function validates inputs and creates a `MultiAssayExperiment` object.
#' If omics data is provided, it converts matrices into `SummarizedExperiment`
#' objects with proper sample alignment. If omics is `NULL`, the function
#' creates an exposure-only object suitable for epidemiological analyses
#' using `run_association()` with `source = "exposures"`.
#'
#' @return A `MultiAssayExperiment` object containing the formatted exposure
#' and optionally omics datasets.
#'
#' @examples
#' # Epi user workflow
#' # so no omics data
#' epi_data <- data.frame(
#'     pm25 = rnorm(10),
#'     outcome = rbinom(10, 1, 0.5),
#'     age = rnorm(10, 45, 10),
#'     row.names = paste0("subj_", 1:10)
#' )
#'
#' codebook <- data.frame(
#'     variable = c("pm25", "outcome", "age"),
#'     category = c("exposure", "outcome", "covariate")
#' )
#'
#' mae <- create_exposomicset(
#'     codebook = codebook,
#'     exposure = epi_data
#' )
#'
#' # Multi-omics workflow
#' tmp <- make_example_data(n_samples = 10)
#'
#' mae <- create_exposomicset(
#'     codebook = tmp$codebook,
#'     exposure = tmp$exposure,
#'     omics = tmp$omics,
#'     row_data = tmp$row_data
#' )
#'
#' @export
create_exposomicset <- function(
  codebook,
  exposure,
  omics = NULL,
  row_data = NULL
) {
    # Validate exposure input
    if (!is.data.frame(exposure)) {
        stop("The 'exposure' argument must be a data frame.")
    }
    if (is.null(rownames(exposure))) {
        stop("Exposure data must have rownames.")
    }

   col_data <- S4Vectors::DataFrame(exposure, row.names = rownames(exposure))

    # Handle epi-only case (no omics)
   if (is.null(omics)) {
     message("No omics data provided. Creating exposure-only exposomicset.")

     # Create a minimal placeholder experiment with 0 features
     placeholder <- SummarizedExperiment::SummarizedExperiment(
       assays = S4Vectors::SimpleList(
         placeholder = matrix(
           nrow = 0,
           ncol = nrow(exposure),
           dimnames = list(NULL, rownames(exposure))
         )
       )
     )

     mae <- MultiAssayExperiment::MultiAssayExperiment(
       experiments = list(.exposures = placeholder),
       colData = col_data,
       metadata = list(codebook = codebook)
     )
     message("MultiAssayExperiment created successfully.")
     return(mae)
   }

    # Validate omics-related inputs
    if (!is.null(row_data) && !is.list(row_data)) {
        stop("The 'row_data' argument must be a list if provided.")
    }

    # Convert a single matrix input into a named list
    if (is.matrix(omics)) {
        omics <- list(single_omic = omics)
    } else if (!is.list(omics)) {
        stop("The 'omics' argument must be a list or a single matrix.")
    }

    if (!is.null(row_data) && length(row_data) != length(omics)) {
        stop("Length of 'row_data' must match 'omics'.")
    }

    # Ensure all datasets in omics are matrices with column names
    message("Ensuring all omics datasets are matrices with column names.")
    omics <- lapply(omics, function(x) {
        if (!is.matrix(x)) x <- as.matrix(x)
        if (is.null(colnames(x))) {
            colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
        }
        x
    })

    # Create row_data if not provided and ensure feature order matches assay rows
    row_data <- row_data %||% lapply(omics, function(df) {
        S4Vectors::DataFrame(row.names = rownames(df))
    })

    # Ensure rowData order matches assay rownames
    row_data <- mapply(function(row_meta, data) {
        row_meta <- row_meta[
            match(rownames(data), rownames(row_meta)), ,
            drop = FALSE
        ]
        row_meta
    }, row_meta = row_data, data = omics, SIMPLIFY = FALSE)

    # Create SummarizedExperiment objects with ordered samples and rowData
    message("Creating SummarizedExperiment objects.")
    experiments <- mapply(
        function(data, row_meta) {
            sample_order <- sort(colnames(data))
            data <- data[, sample_order, drop = FALSE]
            SummarizedExperiment::SummarizedExperiment(
                assays = S4Vectors::SimpleList(counts = data),
                rowData = row_meta
            )
        },
        data = omics,
        row_meta = row_data,
        SIMPLIFY = FALSE
    )

    # Create MultiAssayExperiment
    message("Creating MultiAssayExperiment object.")
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiments,
        colData = col_data,
        metadata = list(codebook = codebook)
    )

    message("MultiAssayExperiment created successfully.")
    return(mae)
}
