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
    resampling = TRUE,  
    resampling_iterations = 50,
    cross_validation = FALSE,  
    k_folds = 5,
    score_thresh=NULL,
    score_quantile=0.1,
    action="add"
) {
  require(tidybulk)
  require(MultiAssayExperiment)
  require(tidyverse)
  
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
            for (exp_name in names(experiments(expomicset))) {
              message("Processing experiment: ", exp_name)
              
              exp <- .update_assay_colData(expomicset, exp_name)
              
              # Skip if too few features
              features_to_test <- exp |> 
                identify_abundant(minimum_counts = min_counts,
                                  minimum_proportion = min_prop) |>
                elementMetadata() |> 
                as.data.frame() |>  
                filter(.abundant == TRUE) |> 
                nrow()
              
              if (features_to_test < 2) {
                warning("Skipping assay ", exp_name, " due to insufficient features.")
                next
              }
              
              # If DESeq2 is used, ensure integer values
              if (method == "DESeq2" && !all(assay(exp) == floor(assay(exp)))) {
                message("Detected non-integer values for DESeq2 in ", exp_name, ". Rounding to nearest integer...")
                assay(exp, abundance_col) <- round(assay(exp, abundance_col), 0)
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
                res <- res |> mutate(model = model_name, exp_name = exp_name)
                sensitivity_df <- bind_rows(sensitivity_df, res)
              }
            }
          }
        }
      }
    }
  }

  
  # Determine stable features 
  feature_stability_df <- sensitivity_df |> 
    .calculate_feature_stability()
  
  
  if(is.null(score_thresh)){
    score_thresh <- quantile(
      feature_stability_df$stability_score,
      score_quantile)
  }else{
    score_thresh <- score_thresh
  }
  
  sum <- feature_stability_df |> 
    group_by(exp_name) |>
    dplyr::reframe(
      n_above=sum(stability_score>score_thresh),
      n=n())
  
  message("Number of features above threshold of ", score_thresh, ":")
  message("-----------------------------------------")
  for(exp_name in unique(feature_stability_df$exp_name)){
    n_above  <- sum |> 
      filter(exp_name==!!exp_name) |>
      pull(n_above);
    n <- sum |> 
      filter(exp_name==!!exp_name) |>
      pull(n);
    message(exp_name, ": ", n_above,"/", n)
  }
  
  message("Sensitivity analysis completed.")
  
  if(action =="add"){
    # Store results in metadata
    metadata(expomicset)$sensitivity_analysis <- list(
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
