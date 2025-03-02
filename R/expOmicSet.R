expOmicSet <- function(var_info, exposure, omics, row_data = NULL) {
  require(MultiAssayExperiment)
  
  # Validate inputs
  if (!is.data.frame(exposure)) stop("The 'exposure' argument must be a data frame.")
  if (!is.null(row_data) && !is.list(row_data)) stop("The 'row_data' argument must be a list if provided.")
  
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
  message("Ensuring all omics datasets are matrices with column names...")
  omics <- lapply(omics, function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
    }
    x
  })
  
  # Create row_data if not provided and ensure feature order matches assay rows
  row_data <- row_data %||% lapply(omics, function(df) {
    DataFrame(row.names = rownames(df))
  })
  
  # Ensure rowData order matches assay rownames
  row_data <- mapply(function(row_meta, data) {
    row_meta <- row_meta[match(rownames(data), rownames(row_meta)), , drop = FALSE]
    row_meta
  }, row_meta = row_data, data = omics, SIMPLIFY = FALSE)
  
  # Create SummarizedExperiment objects with ordered samples and rowData
  message("Creating SummarizedExperiment objects with ordered samples and matching rowData...")
  experiments <- mapply(
    function(data, row_meta) {
      sample_order <- sort(colnames(data))  # Enforce ordered sample names
      data <- data[, sample_order, drop = FALSE]  # Reorder columns
      SummarizedExperiment(
        assays = SimpleList(counts = data),
        rowData = row_meta
      )
    },
    data = omics,
    row_meta = row_data,
    SIMPLIFY = FALSE
  )
  
  # Create MultiAssayExperiment without enforcing sample consistency
  message("Creating MultiAssayExperiment object...")
  mae <- MultiAssayExperiment(
    experiments = experiments,
    colData = DataFrame(exposure),  # Keep all exposure samples without filtering
    metadata = list(var_info = var_info)
  )
  
  message("MultiAssayExperiment created successfully with ordered samples and rowData.")
  return(mae)
}

# --- Second Try ---------
# expOmicSet <- function(var_info, exposure, omics, row_data = NULL) {
#   require(MultiAssayExperiment)
#   
#   # Validate inputs
#   if (!is.data.frame(exposure)) stop("The 'exposure' argument must be a data frame.")
#   if (!is.null(row_data) && !is.list(row_data)) stop("The 'row_data' argument must be a list if provided.")
#   
#   # Convert a single matrix input into a named list
#   if (is.matrix(omics)) {
#     omics <- list(single_omic = omics)
#   } else if (!is.list(omics)) {
#     stop("The 'omics' argument must be a list or a single matrix.")
#   }
#   
#   if (!is.null(row_data) && length(row_data) != length(omics)) {
#     stop("Length of 'row_data' must match 'omics'.")
#   }
#   
#   # Ensure all datasets in omics are matrices with column names
#   message("Ensuring all omics datasets are matrices with column names...")
#   omics <- lapply(omics, function(x) {
#     if (!is.matrix(x)) x <- as.matrix(x)
#     if (is.null(colnames(x))) {
#       colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
#     }
#     x
#   })
#   
#   # Create row_data if not provided
#   row_data <- row_data %||% lapply(omics, function(df) {
#     DataFrame(row.names = rownames(df))
#   })
#   
#   # Create SummarizedExperiment objects with per-omics sample ordering
#   message("Creating SummarizedExperiment objects with ordered samples...")
#   experiments <- mapply(
#     function(data, row_meta) {
#       sample_order <- sort(colnames(data))  # Enforce ordered sample names
#       data <- data[, sample_order, drop = FALSE]  # Reorder columns
#       SummarizedExperiment(
#         assays = SimpleList(counts = data),
#         rowData = row_meta
#       )
#     },
#     data = omics,
#     row_meta = row_data,
#     SIMPLIFY = FALSE
#   )
#   
#   # Create MultiAssayExperiment without enforcing sample consistency
#   message("Creating MultiAssayExperiment object...")
#   mae <- MultiAssayExperiment(
#     experiments = experiments,
#     colData = DataFrame(exposure),  # Keep all exposure samples without filtering
#     metadata = list(var_info = var_info)
#   )
#   
#   message("MultiAssayExperiment created successfully with ordered samples per assay.")
#   return(mae)
# }

# --- First Try ------------

# expOmicSet <- function(var_info, exposure, omics, row_data = NULL) {
#   require(MultiAssayExperiment)
#   
#   # validate inputs
#   if (!is.data.frame(exposure)) stop("The 'exposure' argument must be a data frame.")
#   if (!is.list(omics)) stop("The 'omics' argument must be a list.")
#   if (!is.null(row_data) && !is.list(row_data)) stop("The 'row_data' argument must be a list if provided.")
#   if (!is.null(row_data) && length(row_data) != length(omics)) stop("Length of 'row_data' must match 'omics'.")
#   
#   # ensure all datasets in omics are matrices and have column names
#   message("Ensuring all omics datasets are matrices with column names...")
#   omics <- lapply(omics, function(x) {
#     if (!is.matrix(x)) x <- as.matrix(x)
#     if (is.null(colnames(x))) {
#       colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
#     }
#     x
#   })
#   
#   # create row_data if not provided
#   row_data <- row_data %||% lapply(omics, function(df) {
#     DataFrame(row.names = rownames(df))
#   })
#   
#   # create summarized experiments for each omics dataset
#   message("Creating SummarizedExperiment objects...")
#   experiments <- mapply(
#     function(data, row_meta) {
#       SummarizedExperiment(
#         assays = SimpleList(counts=data),
#         rowData = row_meta
#       )
#     },
#     data = omics,
#     row_meta = row_data,
#     SIMPLIFY = FALSE
#   )
#   
#   # create the MultiAssayExperiment object
#   message("Creating MultiAssayExperiment object with all samples...")
#   mae <- MultiAssayExperiment(
#     experiments = experiments,
#     colData = DataFrame(exposure),  # Include all exposure samples
#     metadata = list(var_info = var_info)
#   )
#   
#   message("MultiAssayExperiment created successfully.")
#   return(mae)
# }
