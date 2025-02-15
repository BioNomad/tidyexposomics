transform_exposure <- function(expOmicSet, transform_method = "best") {
  require(tidyverse)
  require(broom)
  require(ggpubr)
  require(ggsci)
  
  # Extract variables to transform from normality metadata
  variables_to_transform <- metadata(expOmicSet)$normality$norm_df$exposure
  
  if (is.null(variables_to_transform) || length(variables_to_transform) == 0) {
    stop("Please run `check_normality` before transforming the exposure data.")
  }
  
  # Helper function to apply transformations
  apply_transformation <- function(data, columns, method) {
    if (method == "none") {
      data
    } else if (method == "log2") {
      data |> mutate(across(all_of(columns), ~ log2(. + abs(min(.)) + 1)))
    } else if (method == "x_1_3") {
      data |> mutate(across(all_of(columns), ~ (. + abs(min(.)))^(1/3)))
    } else if (method == "sqrt") {
      data |> mutate(across(all_of(columns), ~ sqrt(. + abs(min(.)))))
    } else {
      stop("Unsupported transformation method: ", method)
    }
  }
  
  # Perform transformations based on the method
  if (transform_method == "best") {
    message("Evaluating the best transformation method including no transformation.")
    
    # Separate numeric and non-numeric columns
    col_data <- colData(expOmicSet) |> as.data.frame()
    numeric_data <- col_data |> dplyr::select(all_of(variables_to_transform))
    non_numeric_data <- col_data |> dplyr::select(-all_of(variables_to_transform))
    
    # Apply transformations to numeric columns
    transformations <- list(
      none_trans = apply_transformation(numeric_data, variables_to_transform, "none"),
      log_trans = apply_transformation(numeric_data, variables_to_transform, "log2"),
      x_1_3_trans = apply_transformation(numeric_data, variables_to_transform, "x_1_3"),
      sqrt_trans = apply_transformation(numeric_data, variables_to_transform, "sqrt")
    )
    
    # Evaluate normality for each transformation
    norm_results <- lapply(names(transformations), function(name) {
      transformed <- transformations[[name]]
      transformed |>
        apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
        (\(x) do.call(rbind, x))() |>
        mutate(transformation = name, exposure = colnames(transformed))
    }) |>
      bind_rows()
    
    # Summarize normality results
    norm_summary <- norm_results |>
      group_by(transformation) |>
      dplyr::summarise(
        n = dplyr::n(),
        normal = sum(p.value > 0.05),
        not_normal = sum(p.value <= 0.05),
        percent_normal = normal / n
      ) |>
      arrange(desc(percent_normal))
    
    # Select the best transformation
    best_transformation <- norm_summary$transformation[1]
    if (best_transformation == "none_trans") {
      message("No transformation applied.")
    } else {
      message("Using the ", best_transformation, " transformation.")
    }
    
    # Update colData with the best transformation
    transformed_numeric <- transformations[[best_transformation]]
    updated_col_data <- bind_cols(non_numeric_data, transformed_numeric)
    colData(expOmicSet) <- DataFrame(updated_col_data)
    
    # only save the results for the chosen transformation
    norm_results <- norm_results |> filter(transformation == best_transformation)
    
    # Save results
    metadata(expOmicSet)$transformation <- list(
      norm_df = norm_results,
      norm_summary = norm_summary
    )
    
  } else {
    # Apply a specific transformation
    message("Applying the ", transform_method, " transformation.")
    
    # Separate numeric and non-numeric columns
    col_data <- colData(expOmicSet) |> as.data.frame()
    numeric_data <- col_data |> dplyr::select(all_of(variables_to_transform))
    non_numeric_data <- col_data |> dplyr::select(-all_of(variables_to_transform))
    
    # Transform the numeric data
    transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)
    
    # Combine with non-numeric columns
    updated_col_data <- bind_cols(non_numeric_data, transformed_numeric)
    colData(expOmicSet) <- DataFrame(updated_col_data)
    
    # Perform Shapiro-Wilk test
    norm_results <- transformed_numeric |>
      apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
      (\(x) do.call(rbind, x))() |>
      mutate(exposure = colnames(transformed_numeric))
    
    # Save results
    metadata(expOmicSet)$transformation <- list(
      norm_df = norm_results
    )
  }
  
  return(expOmicSet)
}
