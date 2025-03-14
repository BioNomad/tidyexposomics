#' Associate Latent Factors with an Outcome
#'
#' Tests associations between integration-derived latent factors and an outcome 
#' using correlation or Kruskal-Wallis tests, depending on whether the outcome is categorical.
#' Applies multiple testing correction and optionally stores results in `expomicset` metadata.
#'
#' @param expomicset A \code{MultiAssayExperiment} object containing integration results.
#' @param outcome_var A character string specifying the outcome variable name.
#' @param categorical A logical value indicating whether the outcome is categorical (\code{TRUE}) or continuous (\code{FALSE}). Default is \code{FALSE}.
#' @param p_thresh A numeric value specifying the p-value threshold for significance filtering. Default is \code{0.05}.
#' @param action A character string indicating whether to return results (\code{"get"}) or add them to metadata (\code{"add"}). Default is \code{"add"}.
#'
#' @details
#' The function extracts integration-derived factors from \code{metadata(expomicset)} based on the integration method 
#' (\code{"MOFA"} or \code{"MCIA"}). It then matches factor sample identifiers with those of the outcome and performs 
#' either Spearman correlation (\code{cor.test()}) for continuous outcomes or a Kruskal-Wallis test (\code{kruskal.test()}) 
#' for categorical outcomes. Significant results (\code{p.value < p_thresh}) are retained.
#'
#' @return If \code{action = "add"}, returns the modified \code{expomicset} with significant factors stored in metadata.
#' If \code{action = "get"}, returns a data frame containing:
#' \item{factor}{The latent factor being tested.}
#' \item{outcome}{The outcome variable.}
#' \item{estimate}{The test statistic (correlation coefficient or Kruskal-Wallis chi-squared).}
#' \item{p.value}{The p-value of the test.}
#'
#' @examples
#' \dontrun{
#' results <- associate_factor_outcome(
#'   expomicset = my_expomicset,
#'   outcome_var = "fev_height",
#'   categorical = TRUE,
#'   p_thresh = 0.05,
#'   action = "get"
#' )
#' }
#'
#' @export
associate_factor_outcome <- function(
    expomicset, 
    outcome_var, 
    categorical = FALSE, 
    p_thresh = 0.05,
    action = "add") {
  require(Hmisc)
  require(MultiAssayExperiment)
  require(stats)
  require(broom)
  require(dplyr)
  require(purrr)
  require(tibble)
  
  message("Extracting factors from integration results...")
  
  # Get integration results
  integration_results <- MultiAssayExperiment::metadata(expomicset)$integration_results
  method_used <- integration_results$method
  
  # Extract factors based on integration method
  if (method_used == "MOFA") {
    message("Detected MOFA+ results, extracting factors...")
    factors <- MOFA2::get_factors(integration_results$result)[[1]]
  } else if (method_used == "MCIA") {
    message("Detected MCIA results, extracting global scores...")
    factors <- integration_results$result@global_scores
  } else {
    stop("Unsupported integration method: ", method_used)
  }
  
  # Check if factors are available
  if (is.null(factors)) {
    stop("No factors found in integration results.")
  }
  
  # Extract outcome variable from colData
  outcome_data <- MultiAssayExperiment::colData(expomicset) |> 
    as.data.frame()
  outcome_data <- outcome_data[, outcome_var, drop = FALSE]
  
  # Ensure common samples
  common_samples <- intersect(rownames(factors), rownames(outcome_data))
  if (length(common_samples) < 2) {
    stop("Not enough overlapping samples between factors and outcome for correlation.")
  }
  
  # Subset to common samples
  factors <- factors |> 
    as.data.frame() |> 
    (\(df) df |> 
       dplyr::filter(rownames(df) %in% common_samples))()
  outcome_data <- outcome_data |> 
    as.data.frame() |> 
    (\(df) df |> 
       dplyr::filter(rownames(df) %in% common_samples))()
  
  # Merge factors and outcome into a single data frame
  merged <- factors |> 
    as.data.frame() |> 
    tibble::rownames_to_column("id") |> 
    dplyr::inner_join(outcome_data |> 
                 as.data.frame() |> 
                 tibble::rownames_to_column("id"),
               by = "id")
  
  message("Running correlation or Kruskal-Wallis test...")
  
  # Compute correlation/Kruskal-Wallis test using map2()
  results <- purrr::map2(
    colnames(factors), 
    outcome_var, 
    \(factor, outcome) {
    if (categorical) {
      kruskal.test(merged[[factor]] ~ merged[[outcome]]) |> 
        broom::tidy() |> 
        dplyr::mutate(factor = factor, 
                      outcome = outcome)
    } else {
      cor.test(merged[[factor]], 
               merged[[outcome]],
               method = "spearman",
               exact = FALSE) |> 
        broom::tidy() |>
        dplyr::mutate(factor = factor, 
                      outcome = outcome)
    }
  }) |> 
    dplyr::bind_rows()
  
  # Adjust for multiple testing (Benjamini-Hochberg FDR correction)
  #results <- results |> mutate(FDR = p.adjust(p.value, method = "fdr"))
  
  # Keep only significant results
  results <- results |> 
    dplyr::filter(p.value < p_thresh)
  
  message("Factors with significant association with outcome:")
  for (i in 1:nrow(results)) {
    message(paste0("Factor: ",
                   results$factor[i],
                   ", p-value: ",
                   round(results$p.value[i],digits = 3)))
  }
  if (action=="add"){
    # Store significant factors in metadata
    MultiAssayExperiment::metadata(expomicset)$significant_factors <- results
    
    return(expomicset)
  }else if (action =="get"){
    return(results)
  }else{
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
}
