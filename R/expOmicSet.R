# expOmicSet <- function(var_info, exposure, omics, row_data = NULL) {
#   require(MultiAssayExperiment)
#   # Validate inputs
#   if (!is.data.frame(exposure)) stop("The 'exposure' argument must be a data frame.")
#   if (!is.list(omics)) stop("The 'omics' argument must be a list.")
#   if (!is.null(row_data) && !is.list(row_data)) stop("The 'row_data' argument must be a list if provided.")
#   if (!is.null(row_data) && length(row_data) != length(omics)) stop("Length of 'row_data' must match 'omics'.")
#   
#   # Check common samples
#   common_samples <- rownames(exposure)
#   if (is.null(common_samples)) stop("Row names are missing in the exposure data frame.")
#   
#   for (omics_name in names(omics)) {
#     common_samples <- intersect(common_samples, colnames(omics[[omics_name]]))
#   }
#   if (length(common_samples) == 0) stop("No common samples between exposure and omics data.")
#   
#   # Subset exposure and omics by common samples
#   exposure <- exposure[common_samples, , drop = FALSE]
#   omics <- lapply(omics, function(df) {
#     df <- df[, common_samples, drop = FALSE]
#     as.matrix(df)  # Ensure omics data is a numeric matrix
#   })
#   
#   # Create row_data if not provided
#   row_data <- row_data %||% lapply(omics, function(df) {
#     DataFrame(row.names = rownames(df))
#   })
#   
#   # Create SummarizedExperiment objects for each omics dataset
#   experiments <- mapply(
#     function(data, row_meta) {
#       SummarizedExperiment(
#         assays = SimpleList(data),
#         colData = DataFrame(exposure),
#         rowData = row_meta
#       )
#     },
#     data = omics,
#     row_meta = row_data,
#     SIMPLIFY = FALSE
#   )
# 
#   
#   # Create MultiAssayExperiment
#   mae <- MultiAssayExperiment(
#     experiments = experiments,
#     colData = DataFrame(exposure),  # Central sample metadata
#     metadata = list(var_info=var_info)  # Variable information
#   )
#   
#   
#   return(mae)
# }

expOmicSet <- function(var_info, exposure, omics, row_data = NULL) {
  require(MultiAssayExperiment)
  
  # validate inputs
  if (!is.data.frame(exposure)) stop("The 'exposure' argument must be a data frame.")
  if (!is.list(omics)) stop("The 'omics' argument must be a list.")
  if (!is.null(row_data) && !is.list(row_data)) stop("The 'row_data' argument must be a list if provided.")
  if (!is.null(row_data) && length(row_data) != length(omics)) stop("Length of 'row_data' must match 'omics'.")
  
  # ensure all datasets in omics are matrices and have column names
  message("Ensuring all omics datasets are matrices with column names...")
  omics <- lapply(omics, function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
    }
    x
  })
  
  # create row_data if not provided
  row_data <- row_data %||% lapply(omics, function(df) {
    DataFrame(row.names = rownames(df))
  })
  
  # create summarized experiments for each omics dataset
  message("Creating SummarizedExperiment objects...")
  experiments <- mapply(
    function(data, row_meta) {
      SummarizedExperiment(
        assays = SimpleList(counts=data),
        rowData = row_meta
      )
    },
    data = omics,
    row_meta = row_data,
    SIMPLIFY = FALSE
  )
  
  # create the MultiAssayExperiment object
  message("Creating MultiAssayExperiment object with all samples...")
  mae <- MultiAssayExperiment(
    experiments = experiments,
    colData = DataFrame(exposure),  # Include all exposure samples
    metadata = list(var_info = var_info)
  )
  
  message("MultiAssayExperiment created successfully.")
  return(mae)
}
