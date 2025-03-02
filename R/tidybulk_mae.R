
pivot_sample <- function(x, ...) {
  # Load required libraries
  require(MultiAssayExperiment)
  require(tidyverse)

  if (inherits(x, "MultiAssayExperiment")) {
    # Call your custom function
    colData(x) |>
      as.data.frame() |>
      rownames_to_column(".sample") |>
      as_tibble()
  } else if (inherits(x, "SummarizedExperiment")) {
    # Call the pivot_sample() from the other package
    tidybulk::pivot_sample(x, ...)
  } else {
    stop("Error: pivot_sample() only supports MultiAssayExperiment and SummarizedExperiment objects.")
  }
}
# 
# 
# # Test with MultiAssayExperiment and SummarizedExperiment objects
# expom_enrich |> pivot_sample()
# .update_assay_colData(expom_qc,"Serum Proteomics") |> pivot_sample()
# 
# 
pivot_feature <- function(expOmicSet){
  # Load required libraries
  require(MultiAssayExperiment)
  require(tidyverse)
  require(tidybulk)

  res <- lapply(names(experiments(expOmicSet)),function(exp_name){
    exp <- .update_assay_colData(expOmicSet,exp_name) |>
      pivot_transcript()
  })

  names(res) <- names(experiments(expOmicSet))

  res <- res |>
    bind_rows(.id=".exp_name")
  return(res)
}
# 
# expom_enrich |> pivot_feature()
# 
# detect_target <- function(.data, ...) {
#   # Extract filtering condition as an expression
#   args_expr <- substitute(list(...))[-1]
#   
#   # Extract just column names (ignoring the expression part)
#   args <- sapply(args_expr, function(x) as.character(x[[2]]), USE.NAMES = FALSE)
#   
#   # Convert colData() and rowData() to tibbles for proper column name detection
#   col_data_tbl <- as_tibble(colData(.data), rownames = ".sample")
#   col_data_cols <- colnames(col_data_tbl)
#   
#   row_data_tbls <- lapply(experiments(.data), function(exp) as_tibble(rowData(exp), rownames = ".feature"))
#   row_data_cols <- unique(unlist(lapply(row_data_tbls, colnames)))
#   
#   # Determine where the column exists
#   if (all(args %in% col_data_cols)) {
#     return("colData")
#   } else if (all(args %in% row_data_cols)) {
#     return("rowData")
#   } else {
#     stop("Error: Specified columns do not exist in colData or rowData.")
#   }
# }
# 
# 
# filter.MultiAssayExperiment <- function(.data, ...) {
#   target <- detect_target(.data, ...)
#   
#   if (target == "colData") {
#     # Convert rownames into `.sample` before filtering
#     coldata_tbl <- as_tibble(colData(.data), rownames = ".sample")
#     
#     # Apply filtering
#     coldata_filtered <- coldata_tbl |> filter(...)
#     
#     # Match samples based on `.sample` column
#     matching_samples <- coldata_filtered$.sample
#     .data <- .data[, matching_samples, drop = FALSE]
#     
#   } else if (target == "rowData") {
#     # Subset rowData() across all experiments
#     new_experiments <- lapply(experiments(.data), function(exp) {
#       rowdata_tbl <- as_tibble(rowData(exp), rownames = ".feature")
#       
#       # Apply filtering
#       rowdata_filtered <- rowdata_tbl |> filter(...)
#       matching_features <- rowdata_filtered$.feature
#       
#       # Subset experiment based on features
#       exp[matching_features, , drop = FALSE]
#     })
#     
#     # Recreate the MultiAssayExperiment with filtered experiments
#     .data <- MultiAssayExperiment(
#       experiments = new_experiments,
#       colData = colData(.data),
#       sampleMap = sampleMap(.data),
#       metadata = metadata(.data)
#     )
#   }
#   
#   return(.data)
# }
# 
# 
# 
# 
# a=expom_enrich 
# 
# a |> filter(.sample %in% "s1724")  # Filters based on sample ID
# a |> filter(race %in% "non_black")  # Filters based on colData()
# a |> filter(.feature == "TP53") # Filters based on rowData()

# ---------- Test where I use lists -----------



# tidybulk_mae <- function(expOmicSet) {
#   library(MultiAssayExperiment)
#   library(tidybulk)
#   library(dbplyr)
#   library(dplyr)
#   res <- lapply(names(experiments(expOmicSet)),
#                 function(name){
#                   res=.update_assay_colData(a,name) |> 
#                     tidybulk();
#                   gc();
#                   return(res)})
# }


tidybulk_mae <- function(expOmicSet) {
  library(MultiAssayExperiment)
  library(tidybulk)
  library(dplyr)
  library(dbplyr)
  
  results_env <- new.env()  # Store references in an environment
  
  for (name in names(experiments(expOmicSet))) {
    message("Processing: ", name)
    
    # Process each experiment
    res <- .update_assay_colData(expOmicSet, name) |> tidybulk()
    
    # Store reference
    results_env[[name]] <- res
    
    # Remove intermediate object & free memory
    rm(res)
    gc()
  }
  
  return(results_env)  # Return environment with references
}

#x=a |> tidybulk_mae()
