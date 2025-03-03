adjust_assays <- function(
    expomicset, 
    outcome, 
    covariates, 
    minimum_counts = 0,
    minimum_proportion = 0.7,
    scaling_method = "none",
    skip_identify_abundant = NULL
) {
  require(tidybulk)
  require(MultiAssayExperiment)
  require(tidyverse)
  
  message("Adjusting assays...")
  adjusted_experiments <- list()
  
  for (assay_name in names(experiments(expomicset))) {
    message("Processing assay: ", assay_name)
    
    # Update colData for the assay
    assay <- .update_assay_colData(expomicset, assay_name)
    
    # Create batch variable in colData
    coldata <- as.data.frame(colData(assay))
    coldata <- coldata |>
      mutate(batch = pmap_chr(across(all_of(covariates)), ~ paste(..., sep = "_")))
    
    # Update colData with batch
    colData(assay) <- DataFrame(coldata)
    
    # Apply `identify_abundant` only if not skipped
    if (!is.null(skip_identify_abundant) && assay_name %in% skip_identify_abundant) {
      message("Skipping `identify_abundant` for assay: ", assay_name)
      assay <- assay |>
        identify_abundant(
          minimum_counts = 0,
          minimum_proportion = 0)
    } else {
      message("Applying `identify_abundant` for assay: ", assay_name)
      assay <- assay |>
        identify_abundant(minimum_counts = minimum_proportion,
                          minimum_proportion = minimum_proportion)
    }
    
    # Always apply `scale_abundance`
    message("Applying `scale_abundance` for assay: ", assay_name)
    assay <- assay |>
      scale_abundance(method = scaling_method)
    
    # Adjust abundance
    message("Adjusting abundance for assay: ", assay_name)
    assay <- assay |>
      adjust_abundance(as.formula(paste("~", paste(outcome, " + batch", sep = ""))), 
                       method = "combat")
    
    # Remove colData from assay to keep expomicset central colData as the primary source
    colData(assay) <- NULL
    
    # Save adjusted assay
    adjusted_experiments[[assay_name]] <- assay
  }
  
  # Update experiments in the expomicset
  experiments(expomicset) <- ExperimentList(adjusted_experiments)
  
  message("Adjusted variation due to: ", paste(covariates, collapse = ", "))
  return(expomicset)
}


# adjust_assays <- function(expomicset, outcome, covariates) {
#   require(tidybulk)
#   require(MultiAssayExperiment)
#   require(tidyverse)
#   
#   # Initialize a list to store adjusted experiments
#   adjusted_experiments <- list()
#   
#   # Iterate through all assays in the expomicset
#   message("Adjusting assays...")
#   for (assay_name in names(experiments(expomicset))) {
#     message("Processing assay: ", assay_name)
#     
#     assay <- experiments(expomicset)[[assay_name]]
#     
#     # Extract colData for the assay's samples
#     assay_samples <- colnames(assay)
#     global_coldata <- as.data.frame(colData(expomicset))
#     coldata <- global_coldata[rownames(global_coldata) %in% assay_samples, , drop = FALSE]
#     
#     # Ensure the sample order matches
#     coldata <- coldata[match(assay_samples, rownames(coldata)), , drop = FALSE]
#     
#     # Add a check to ensure the order is correct
#     if (!identical(rownames(coldata), assay_samples)) {
#       stop("Sample order mismatch detected in assay: ", assay_name, 
#            "\nEnsure the samples in colData are aligned with the assay samples.")
#     }
#     
#     # Ensure required covariates are present in coldata
#     if (!all(covariates %in% colnames(coldata))) {
#       stop("Missing one or more covariates in colData for assay: ", assay_name)
#     }
#     
#     # Create batch variable
#     message("Creating batch variable for assay: ", assay_name)
#     coldata <- coldata |>
#       mutate(batch = pmap_chr(across(all_of(covariates)), ~ paste(..., sep = "_")))
#     
#     # Update colData in the assay
#     colData(assay) <- DataFrame(coldata)
#     
#     # Adjust abundance
#     message("Adjusting abundance for assay: ", assay_name)
#     assay <- assay |>
#       identify_abundant() |>
#       scale_abundance() |>
#       adjust_abundance(as.formula(paste("~", paste(outcome, " + batch", sep = ""))))
#     
#     # Store adjusted assay in the list
#     adjusted_experiments[[assay_name]] <- assay
#   }
#   
#   # Update experiments in the expomicset
#   experiments(expomicset) <- ExperimentList(adjusted_experiments)
#   
#   message("Adjusted variation due to: ", paste(covariates, collapse = ", "))
#   return(expomicset)
# }



# require(tidybulk)
# 
# a <-  c("female", "income5", "black")
# 
# b <- experiments(expom_10)[["cd4_rna"]]
# 
# c <- colData(b) |> 
#   as.data.frame() |> 
#   mutate(batch = pmap_chr(across(all_of(a)), ~ paste(..., sep = "_")))
# 
# colData(b) <- DataFrame(c)
# 
# d <- b |> 
#   identify_abundant() |> 
#   scale_abundance() |>
#   adjust_abundance(~ pftfev1fvc_actual + batch)



# cd16_rna |> 
#   t() |> 
#   as.data.frame() |> 
#   dplyr::select(UTY) |>
#   rownames_to_column("aw_id") |> 
#   inner_join(aw_dataset_cv2_filt |> 
#                dplyr::select(female) |> 
#                rownames_to_column("aw_id"),
#              by="aw_id") |> 
#   view()
