#' Create an expomicset Object
#'
#' Constructs a `MultiAssayExperiment` object from exposure data and omics datasets,
#' ensuring proper formatting and alignment of samples and features.
#'
#' @param codebook A data frame containing variable information metadata.
#' @param exposure A data frame containing exposure data,
#' with rows as samples and columns as variables.
#' @param omics A list of matrices or a single matrix representing omics data.
#'  Each matrix should have samples as columns and features as rows.
#' @param row_data An optional list of `DataFrame` objects providing
#' feature metadata for each omics dataset. If `NULL`,
#' row metadata is generated automatically. Default is `NULL`.
#'
#' @details
#' The function validates inputs, converts `omics` into a list if necessary,
#' ensures all datasets are matrices with column names,
#' and creates `SummarizedExperiment` objects for each omics dataset.
#' It then constructs a `MultiAssayExperiment` object
#' with exposure data in `colData` and variable information stored in metadata.
#'
#' @return A `MultiAssayExperiment` object containing the formatted exposure
#' and omics datasets.
#'
#' @examples
#'
#' # make the example data
#' tmp <- make_example_data(n_samples = 10)
#'
#' # create the MultiAssayExperiment Object
#' mae <- create_expomicset(
#'   codebook = tmp$codebook,
#'   exposure = tmp$exposure,
#'   omics = tmp$omics,
#'  row_data = tmp$row_data
#' )
#'
#' @export
create_expomicset <- function(codebook,
                              exposure,
                              omics,
                              row_data = NULL) {

  # Validate inputs
  if (!is.data.frame(exposure)) {
    stop("The 'exposure' argument must be a data frame.")
  }
  if (!is.null(row_data) && !is.list(row_data)){
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
      match(rownames(data), rownames(row_meta)), , drop = FALSE]
    row_meta
  }, row_meta = row_data, data = omics, SIMPLIFY = FALSE)

  # Create SummarizedExperiment objects with ordered samples and rowData
  message("Creating SummarizedExperiment objects.")
  experiments <- mapply(
    function(data, row_meta) {
      sample_order <- sort(colnames(data))  # Enforce ordered sample names
      data <- data[, sample_order, drop = FALSE]  # Reorder columns
      SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = data),
        rowData = row_meta
      )
    },
    data = omics,
    row_meta = row_data,
    SIMPLIFY = FALSE
  )

  col_data <- S4Vectors::DataFrame(exposure)
  if (is.null(rownames(col_data))) {
    stop("Exposure data must have rownames.")
  }

  # Create MultiAssayExperiment without enforcing sample consistency
  message("Creating MultiAssayExperiment object.")
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    colData = col_data,  # Keep all exposure samples without filtering
    metadata = list(codebook = codebook)
  )

  message("MultiAssayExperiment created successfully.")
  return(mae)
}

