#' Run exposure-wide association study (ExWAS) across all omics layers
#'
#' This function performs an association analysis between a specified outcome
#' and all features across multiple omics layers in a `MultiAssayExperiment` object.
#' The association is performed using generalized linear models (GLMs),
#' optionally adjusting for covariates.
#'
#' Optionally, the top `n` most variable features can be selected from each omics
#' layer prior to modeling to reduce computation time.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics data and colData.
#' @param outcome A character string naming the outcome variable to model.
#' @param covariates A character vector of covariate names to adjust for (optional).
#' @param family A string specifying the GLM family: `"gaussian"` (default) or `"binomial"`.
#' @param top_n_omics Integer. Number of top variable features to select per omics layer. If `NULL`, all features are used.
#' @param correction_method Method used for p-value adjustment (e.g., `"fdr"`, `"bonferroni"`). Default is `"fdr"`.
#' @param action A string indicating whether to `"add"` results to the metadata (default) or `"get"` them as a list.
#'
#' @return If `action = "add"`, returns the original `MultiAssayExperiment` with results
#' stored in `metadata(expomicset)$exwas_all`. If `action = "get"`, returns a list with:
#' - `results_df`: A data frame of association results
#' - `covariates`: Covariates used in the model
#' - `selected_features`: List of selected features per omics layer
#' - `data`: The combined design matrix used in modeling
#'
#' @details
#' - This function scales numeric variables in `colData()`.
#' - Features are named using the pattern `omicsLayer_feature`, and stripped later for reporting.
#' - Handles both continuous and binary outcomes via `family`.
#'
#' @examples
#' # Run ExWAS on a MultiAssayExperiment object
#' # result <- associate_all_outcome(mae, outcome = "BMI", family = "gaussian")
#'
#' @export

associate_all_outcome <- function(expomicset,
                                  outcome,
                                  covariates = NULL,
                                  family = "gaussian",
                                  top_n_omics = NULL,
                                  correction_method = "fdr",
                                  action="add") {

  # Extract and preprocess colData
  message("Extracting exposure data...")
  data <- MultiAssayExperiment::colData(expomicset)  |>
    as.data.frame() |>
    dplyr::mutate_if(is.numeric, ~ scale(.))

  # Log2-transform omics data
  log_2_omics <- .log2_multiassay(expomicset)

  # Check if top_n_omics is provided
  if(!is.null(top_n_omics)){
    selected_features <- .top_var_multiassay(
      expomicset = log_2_omics,
      n = top_n_omics
    )
  } else{
    # Select all features if top_n_omics is not provided
    selected_features <- lapply(
      experiments(log_2_omics),
      function(x) {
        rownames(x)
      }
    )
  }

  # Scale omics data
  scaled_omics <- .scale_multiassay(
    expomicset = log_2_omics,
    log2 = FALSE
  )

  # Extract scaled omics data
  omics_df <- lapply(names(experiments(scaled_omics)),function(name){
    # Extract omics data
    se=experiments(scaled_omics)[[name]]

    # Filter omics data based on selected features
    se=se[rownames(se) %in% selected_features[[name]],]

    # Convert assay data to data frame
    df=se |>
      assay() |>
      t() |>
      as.data.frame()
    names(df)=
      paste0(name,"_",names(df))
    names(df)=gsub(" |-","_",names(df))

    # Create a unique identifier for each row
    df=df |>
      rownames_to_column("id")
  }) |>
    (\(lst){
      # Name list
      names(lst) <- names(MultiAssayExperiment::experiments(scaled_omics))
      lst
    })() |>
    # Combine all omics data into a single data frame
    purrr::reduce(full_join, by = "id")

  # Merge omics data with exposure data
  data <- data |>
    dplyr::mutate(id = rownames(data)) |>
    dplyr::left_join(omics_df, by = "id")

  if(ncol(data)>5000){
    message("Warning: The number of features is greater than 5000. This may take a while.")
    message("Consider using the top_n_omics parameter to reduce the number of features.")
  }

  # Initialize results list
  results <- list()

  # If family is binomial convert to factor and ensure there are two unique values
  if (family == "binomial") {
    if (!is.factor(data[[outcome]])) {
      data[[outcome]] <- as.factor(data[[outcome]])
    }
    if (length(unique(data[[outcome]])) != 2) {
      stop("Outcome variable must have exactly two unique values for binomial family.")
    }
  }

  # Iterate over each column
  message("Performing ExWAS...")
  for (col in colnames(data)) {
    # Construct formula
    formula <- as.formula(
      paste(outcome, "~",
            col,
            if (!is.null(covariates)) paste("+", paste(covariates, collapse = " +")) else "")
    )

    # Fit the generalized linear model
    model <- tryCatch(
      glm(formula, data = data, family = family),
      error = function(e) {
        message("   * Model fitting failed for ", col, ": ", e$message)
        return(NULL)
      }
    )

    # Skip if model failed to fit
    if (is.null(model)) next

    # Extract model statistics
    model_summary <- broom::tidy(model) |>
      dplyr::filter(term == col)

    # Add column name
    model_summary <- model_summary |>
      dplyr::mutate(var = col)

    # Append to results
    results[[col]] <- model_summary
  }

  # Combine results into a single data frame
  message("Combining results...")
  results_df <- dplyr::bind_rows(results)

  # Adjust p-values for multiple testing
  message("Applying multiple testing correction...")
  results_df <- results_df |>
    dplyr::mutate(p_adjusted = p.adjust(p.value, method = correction_method)) |>
    dplyr::mutate(depend = outcome)

  annot_df <- MultiAssayExperiment::metadata(expomicset)$var_info |>
    dplyr::select(variable, category) |>
    bind_rows(
      selected_features |>
        imap(~ data.frame(
          category = .y,
          variable = .x
        )) |>
        bind_rows() |>
        dplyr::mutate(category = gsub(" ","_",category),
                      variable=paste(category,variable,sep="_"),
                      variable=gsub(" |-","_",variable))
    )

  # Reorder columns for clarity
  results_df <- results_df |>
    dplyr::select(var,
                  depend,
                  estimate,
                  std.error,
                  statistic,
                  p.value,
                  p_adjusted,
                  dplyr::everything()) |>
    dplyr::rename("outcome" = "depend") |>
    dplyr::left_join(annot_df,
                     c("var" = "variable"))
  # dplyr::left_join(MultiAssayExperiment::metadata(expomicset)$var_info,
  #                  by = c("var" = "variable"))

  # Remove omics names from variable names
  names_to_remove <- expomicset |>
    MultiAssayExperiment::experiments() |>
    names() |>
    (\(chr) gsub(" ","_",chr))() |>
    (\(chr) paste0(chr,"_"))() |>
    (\(chr) paste0(chr,collapse = "|"))()

  # Remove omics names from variable names
  results_df <- results_df |>
    dplyr::mutate(var = gsub(names_to_remove,"",var))

  message("ExWAS analysis completed.")

  if(action=="add"){
    # Save results in metadata
    MultiAssayExperiment::metadata(expomicset)$exwas_all <- list(
      results_df=results_df,
      covariates=covariates)

    return(expomicset)
  }else if (action=="get"){
    return(list(
      results_df=results_df,
      covariates=covariates,
      selected_features=selected_features,
      data=data
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
