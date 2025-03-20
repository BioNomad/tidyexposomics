#' Perform Sensitivity Analysis for Differential Abundance
#'
#' Runs differential abundance testing across multiple models, statistical methods,
#' scaling approaches, and filtering criteria to assess feature stability.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param base_formula A formula specifying the base model for differential abundance analysis.
#' @param abundance_col A character string specifying the assay column to use for abundance values. Default is `"counts"`.
#' @param methods A character vector of methods for differential abundance testing (e.g., `"limma_voom"`, `"DESeq2"`, `"edgeR_quasi_likelihood"`).
#' @param scaling_methods A character vector of normalization methods to apply (e.g., `"none"`, `"TMM"`, `"quantile"`).
#' @param min_counts_range A numeric vector of minimum count thresholds to test. Default is `c(1, 5, 10)`.
#' @param min_proportion_range A numeric vector of minimum sample proportion thresholds to test. Default is `c(0.1, 0.2, 0.3)`.
#' @param contrasts A character vector specifying contrasts for the differential analysis. Default is `NULL`.
#' @param covariates_to_remove A character vector of covariates to iteratively remove from the model for sensitivity testing. Default is `NULL`.
#' @param score_thresh A numeric value specifying a fixed threshold for feature stability. If `NULL`, it is determined based on `score_quantile`.
#' @param score_quantile A numeric value specifying the quantile threshold for determining stable features. Default is `0.1`.
#' @param action A character string specifying whether to store (`"add"`) or return (`"get"`) the results. Default is `"add"`.
#'
#' @details
#' This function:
#' - Iterates over multiple **model specifications** by removing covariates.
#' - Tests different **differential abundance methods** (`limma_voom`, `DESeq2`, `edgeR_quasi_likelihood`).
#' - Evaluates **scaling approaches** (`none`, `TMM`, `quantile`).
#' - Filters features based on **minimum count** and **sample proportion** thresholds.
#' - **Calculates feature stability**, measuring how often a feature is identified across conditions.
#' - Determines a **stability threshold** using either `score_thresh` or the `score_quantile`.
#' - **Output Handling**:
#'   - `"add"`: Stores results in `metadata(expomicset)$sensitivity_analysis`.
#'   - `"get"`: Returns a list containing the sensitivity results.
#'
#' @return A `MultiAssayExperiment` object with sensitivity analysis results added to metadata (if `action = "add"`) or a list with:
#' \item{sensitivity_df}{A dataframe of differential abundance results across all tested conditions.}
#' \item{feature_stability}{A dataframe summarizing stability scores for each feature.}
#' \item{score_thresh}{The threshold used for stable feature selection.}
#'
#' @examples
#' \dontrun{
#' expom <- run_sensitivity_analysis(
#'   expomicset = expom,
#'   base_formula = ~ condition + batch,
#'   methods = c("limma_voom", "DESeq2"),
#'   min_counts_range = c(5, 10),
#'   action = "add"
#' )
#'
#' sensitivity_results <- run_sensitivity_analysis(
#'   expomicset = expom,
#'   base_formula = ~ condition + batch,
#'   action = "get"
#' )
#' }
#'
#' @export
run_sensitivity_analysis <- function(
    expomicset,
    base_formula,
    abundance_col = "counts",
    methods = c("limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
    scaling_methods = c("none", "TMM", "quantile"),
    min_counts_range = c(1, 5, 10),
    min_proportion_range = c(0.1, 0.2, 0.3),
    contrasts = NULL,
    covariates_to_remove = NULL,
    pval_col = "adj.P.Val",
    logfc_col = "logFC",
    pval_threshold = 0.05,
    logFC_threshold = log2(1.5),
    score_thresh=NULL,
    score_quantile=0.1,
    action="add"
) {

  message("Running sensitivity analysis for differential abundance...")

  # Extract all terms in the base model
  base_terms <- all.vars(base_formula)

  # Generate all models by removing each covariate one at a time
  model_list <- list()
  model_list[["Full Model"]] <- base_formula

  if (!is.null(covariates_to_remove)) {
    for (covar in covariates_to_remove) {
      reduced_terms <- setdiff(base_terms, covar)
      if (length(reduced_terms) > 1) {
        reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
        model_list[[paste("Without", covar)]] <- reduced_formula
      }
    }
  }

  # Initialize results dataframe
  sensitivity_df <- data.frame()

  for (model_name in names(model_list)) {
    formula <- model_list[[model_name]]

    message("Testing model: ", model_name, " | Formula: ", formula)

    for (method in methods) {
      for (scaling in scaling_methods) {
        for (min_counts in min_counts_range) {
          for (min_prop in min_proportion_range) {

            message("Testing method: ", method,
                    " | Scaling: ", scaling,
                    " | Min Counts: ", min_counts,
                    " | Min Proportion: ", min_prop)

            # Iterate over each experiment in expomicset
            for (exp_name in names(MultiAssayExperiment::experiments(expomicset))) {
              message("Processing experiment: ", exp_name)

              exp <- .update_assay_colData(expomicset, exp_name)

              # Skip if too few features
              features_to_test <- exp |>
                tidybulk::identify_abundant(minimum_counts = min_counts,
                                  minimum_proportion = min_prop) |>
                S4Vectors::elementMetadata() |>
                as.data.frame() |>
                dplyr::filter(.abundant == TRUE) |>
                nrow()

              if (features_to_test < 2) {
                warning("Skipping assay ", exp_name, " due to insufficient features.")
                next
              }

              # If DESeq2 is used, ensure integer values
              if (method == "DESeq2" && !all(assay(exp) == floor(assay(exp)))) {
                message("Detected non-integer values for DESeq2 in ", exp_name, ". Rounding to nearest integer...")
                SummarizedExperiment::assay(exp, abundance_col) <- round(SummarizedExperiment::assay(exp, abundance_col), 0)
              }

              # Call helper function
              res <- .run_se_differential_abundance(
                se = exp,
                formula = formula,
                abundance_col = abundance_col,
                method = method,
                scaling_method = scaling,
                min_counts = min_counts,
                min_proportion = min_prop,
                contrasts = contrasts
              )

              # Append results
              if (!is.null(res)) {
                res <- res |>
                  dplyr::mutate(model = model_name,
                                exp_name = exp_name)
                sensitivity_df <- sensitivity_df |>
                  dplyr::bind_rows(res)
              }
            }
          }
        }
      }
    }
  }


  # Determine stable features
  feature_stability_df <- sensitivity_df |>
    .calculate_feature_stability(pval_col = pval_col,
                                 logfc_col = logfc_col,
                                 pval_threshold = pval_threshold,
                                 logFC_threshold = logFC_threshold)


  if(is.null(score_thresh)){
    score_thresh <- quantile(
      feature_stability_df$stability_score,
      score_quantile)
  }else{
    score_thresh <- score_thresh
  }

  sum <- feature_stability_df |>
    dplyr::group_by(exp_name) |>
    dplyr::reframe(
      n_above=sum(stability_score>score_thresh),
      n=n())

  message("Number of features above threshold of ", score_thresh, ":")
  message("-----------------------------------------")
  for(exp_name in unique(feature_stability_df$exp_name)){
    n_above  <- sum |>
      dplyr::filter(exp_name==!!exp_name) |>
      dplyr::pull(n_above);
    n <- sum |>
      dplyr::filter(exp_name==!!exp_name) |>
      dplyr::pull(n);
    message(exp_name, ": ", n_above,"/", n)
  }

  message("Sensitivity analysis completed.")

  if(action =="add"){
    # Store results in metadata
    MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis <- list(
      sensitivity_df = sensitivity_df,
      feature_stability = feature_stability_df,
      score_thresh = score_thresh
    )
    return(expomicset)
  }else if (action=="get"){
    return(list(
      sensitivity_df = sensitivity_df,
      feature_stability = feature_stability_df,
      score_thresh = score_thresh
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }

}
