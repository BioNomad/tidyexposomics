#' Transform Exposure Data in a MultiAssayExperiment
#'
#' Applies specified or optimal transformations to exposure data in a `MultiAssayExperiment`.
#'
#' @param expomicset A `MultiAssayExperiment` object.
#' @param cols_of_interest Character vector of exposure variables to transform. If `NULL`, all numeric exposures are considered.
#' @param transform_method Transformation method to apply. Options:
#'   - `"none"`: No transformation.
#'   - `"log2"`: Log2 transformation.
#'   - `"x_1_3"`: Cube root transformation.
#'   - `"sqrt"`: Square root transformation.
#'   - `"best"` (default): Selects the transformation maximizing normality based on the Shapiro-Wilk test.
#'
#' @details
#' - If `transform_method = "best"`, the function tests multiple transformations and selects the one that maximizes normality.
#' - Normality is assessed using the Shapiro-Wilk test.
#' - The function updates `colData()` with transformed values.
#' - Transformation results, including normality tests, are stored in `metadata(expomicset)$transformation`.
#'
#' @return A `MultiAssayExperiment` object with transformed exposure data.
#'
#' @examples
#' \dontrun{
#' expom <- transform_exposure(expom, transform_method = "best")
#' expom <- transform_exposure(expom, cols_of_interest = c("Exposure1", "Exposure2"), transform_method = "log2")
#' }
#'
#' @export
transform_exposure <- function(
    expomicset,
    cols_of_interest=NULL,
    transform_method = "best") {

  if(!is.null(cols_of_interest)){
    variables_to_transform <- cols_of_interest
  }else{
    # Extract variables to transform from normality metadata
    variables_to_transform <- MultiAssayExperiment::metadata(expomicset)$normality$norm_df$exposure
  }

  if (is.null(variables_to_transform) || length(variables_to_transform) == 0) {
    stop("Please run `check_normality` before transforming the exposure data.")
  }


  # Helper function to apply transformations
  apply_transformation <- function(data, columns, method) {
    if (method == "none") {
      data
    } else if (method == "log2") {

      # log2 transformation with minimum value adjustment
      data |>
        dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ log2(. + abs(min(.)) + 1)))

    } else if (method == "x_1_3") {

      # x^(1/3) transformation with minimum value adjustment
      data |>
        dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ (. + abs(min(.)))^(1/3)))

    } else if (method == "sqrt") {

      # sqrt transformation with minimum value adjustment
      data |>
        dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ sqrt(. + abs(min(.)))))

    } else {
      stop("Unsupported transformation method: ", method)
    }
  }

  # Confirm that variables are numeric
  col_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame()
  numeric_cols <- col_data |>
    dplyr::select_if(~is.numeric(.)) |>
    colnames()
  variables_to_transform <- variables_to_transform[variables_to_transform %in% numeric_cols]

  # Perform transformations based on the method
  if (transform_method == "best") {
    message("Evaluating the best transformation method including no transformation.")

    # Separate numeric and non-numeric columns
    col_data <- MultiAssayExperiment::colData(expomicset) |>
      as.data.frame()
    numeric_data <- col_data |>
      dplyr::select(dplyr::all_of(variables_to_transform))
    non_numeric_data <- col_data |>
      dplyr::select(-dplyr::all_of(variables_to_transform))

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
        apply(2, function(x) shapiro.test(x) |>
                broom::tidy()) |>
        (\(x) do.call(rbind, x))() |>
        dplyr::mutate(transformation = name,
               exposure = colnames(transformed))
    }) |>
      bind_rows()

    # Summarize normality results
    norm_summary <- norm_results |>
      dplyr::group_by(transformation) |>
      dplyr::summarise(
        n = dplyr::n(),
        normal = sum(p.value > 0.05),
        not_normal = sum(p.value <= 0.05),
        percent_normal = normal / n
      ) |>
      dplyr::arrange(dplyr::desc(percent_normal))

    # Select the best transformation
    best_transformation <- norm_summary$transformation[1]
    if (best_transformation == "none_trans") {
      message("No transformation applied.")
    } else {
      message("Using the ", best_transformation, " transformation.")
    }

    # Update colData with the best transformation
    transformed_numeric <- transformations[[best_transformation]]
    updated_col_data <- non_numeric_data |>
      dplyr::bind_cols(transformed_numeric)
    MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)

    # only save the results for the chosen transformation
    norm_results <- norm_results |>
      dplyr::filter(transformation == best_transformation)

    # Save results
    MultiAssayExperiment::metadata(expomicset)$transformation <- list(
      norm_df = norm_results,
      norm_summary = norm_summary
    )

  } else {
    # Apply a specific transformation
    message("Applying the ", transform_method, " transformation.")

    # Separate numeric and non-numeric columns
    col_data <- MultiAssayExperiment::colData(expomicset) |>
      as.data.frame()
    numeric_data <- col_data |>
      dplyr::select(dplyr::all_of(variables_to_transform))
    non_numeric_data <- col_data |>
      dplyr::select(-dplyr::all_of(variables_to_transform))

    # Transform the numeric data
    transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)

    # Combine with non-numeric columns
    updated_col_data <- non_numeric_data |>
      bind_cols(transformed_numeric)
    MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)

    # Perform Shapiro-Wilk test
    norm_results <- transformed_numeric |>
      dplyr::select_if(~ length(unique(na.omit(.))) >= 3) |>
      apply(2, function(x) shapiro.test(x) |>
              broom::tidy()) |>
      (\(x) do.call(rbind, x))() |>
      dplyr::mutate(exposure = colnames(transformed_numeric))

    # Summarize normality results
    norm_summary <- norm_results |>
      dplyr::summarise(
        "Normal" = sum(p.value > 0.05),
        "Not Normal" = sum(p.value <= 0.05)
      ) |>
      t() |>
      as.data.frame() |>
      rownames_to_column("var") |>
      setNames(c("var","value"))



    # Save results
    MultiAssayExperiment::metadata(expomicset)$transformation <- list(
      norm_df = norm_results,
      norm_summary = norm_summary
    )
  }

  return(expomicset)
}
