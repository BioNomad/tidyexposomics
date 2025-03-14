#' Perform Exposure-Wide Association Study (ExWAS) Analysis
#'
#' Fits generalized linear models (GLMs) to assess associations between multiple exposures 
#' and an outcome, adjusting for covariates. Applies multiple testing correction and 
#' optionally stores results in `expomicset` metadata.
#'
#' @param expomicset An \code{expomicset} object containing exposure and outcome data.
#' @param exposures A character vector of exposure variable names.
#' @param outcome A character string specifying the outcome variable name.
#' @param covariates An optional character vector of covariate names. Default is \code{NULL}.
#' @param family A character string specifying the GLM family (e.g., \code{"gaussian"}, \code{"binomial"}). Default is \code{"gaussian"}.
#' @param correction_method A character string specifying the multiple testing correction method (e.g., \code{"fdr"}, \code{"bonferroni"}). Default is \code{"fdr"}.
#' @param action A character string indicating whether to return results (\code{"get"}) or add them to metadata (\code{"add"}). Default is \code{"add"}.
#'
#' @details
#' The function iterates over exposures, fits a GLM for each, extracts model statistics, 
#' adjusts p-values for multiple testing, and optionally stores results in `expomicset` metadata. 
#' The function handles errors in model fitting gracefully and logs progress with messages.
#'
#' @return If \code{action = "add"}, returns the modified \code{expomicset} with ExWAS results stored in metadata.
#' If \code{action = "get"}, returns a list containing:
#' \item{results_df}{A data frame with model estimates, standard errors, p-values, and adjusted p-values.}
#' \item{covariates}{The covariates used in the models.}
#'
#' @examples
#' \dontrun{
#' results <- associate_exposure_outcome(
#'   expomicset = my_expomicset,
#'   exposures = c("pm25", "lead", "ozone"),
#'   outcome = "lung_function",
#'   covariates = c("age", "sex", "smoking"),
#'   family = "gaussian",
#'   correction_method = "fdr",
#'   action = "get"
#' )
#' }
#'
#' @export
associate_exposure_outcome <- function(expomicset, 
                          exposures, 
                          outcome,
                          covariates = NULL,
                          family = "gaussian",
                          correction_method = "fdr",
                          action="add") {
  
  # Extract and preprocess colData
  message("Extracting exposure data...")
  data <- MultiAssayExperiment::colData(expomicset)  |> 
    as.data.frame() |> 
    dplyr::mutate_if(is.numeric, ~ scale(.))  
  
  # Initialize results list
  results <- list()
  
  # Iterate over each exposure
  message("Performing ExWAS...")
  for (exposure in exposures) {
    # Construct formula
    formula <- as.formula(
      paste(outcome, "~", 
            exposure, 
            if (!is.null(covariates)) paste("+", paste(covariates, collapse = " +")) else "")
    )
    
    # Fit the generalized linear model
    model <- tryCatch(
      glm(formula, data = data, family = family),
      error = function(e) {
        message("   * Model fitting failed for ", exposure, ": ", e$message)
        return(NULL)
      }
    )
    
    # Skip if model failed to fit
    if (is.null(model)) next
    
    # Extract model statistics
    model_summary <- broom::tidy(model) |>
      dplyr::filter(term == exposure)  
    
    # Add exposure name
    model_summary <- model_summary |>
      dplyr::mutate(exposure = exposure)
    
    # Append to results
    results[[exposure]] <- model_summary
  }
  
  # Combine results into a single data frame
  message("Combining results...")
  results_df <- dplyr::bind_rows(results)
  
  # Adjust p-values for multiple testing
  message("Applying multiple testing correction...")
  results_df <- results_df |>
    dplyr::mutate(p_adjusted = p.adjust(p.value, method = correction_method)) |> 
    dplyr::mutate(depend = outcome)
  
  # Reorder columns for clarity
  results_df <- results_df |>
    dplyr::select(exposure,
                  depend, 
                  estimate, 
                  std.error,
                  statistic,
                  p.value, 
                  p_adjusted,
                  dplyr::everything()) |> 
    dplyr::rename("outcome" = "depend") |> 
    dplyr::left_join(MultiAssayExperiment::metadata(expomicset)$var_info,
              by = c("exposure" = "variable"))
  
  message("ExWAS analysis completed.")
  
  if(action=="add"){
    # Save results in metadata
    MultiAssayExperiment::metadata(expomicset)$exwas <- list(
      results_df=results_df,
      covariates=covariates)
    
    return(expomicset)
  }else if (action=="get"){
    return(list(
      results_df=results_df,
      covariates=covariates
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
