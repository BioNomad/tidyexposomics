#' Transform Exposure Variables for Normality
#'
#' Applies a transformation to selected numeric exposure variables in the `colData` of a
#' `MultiAssayExperiment` to improve their normality (e.g., log, Box-Cox, sqrt).
#' Transformation results and normality statistics are stored in metadata for tracking.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure variables in `colData`.
#' @param exposure_cols Optional character vector of exposure variable names to transform.
#'   If `NULL`, uses exposures found in `metadata(expomicset)$quality_control$normality$norm_df$exposure`.
#' @param transform_method Character. Transformation method to apply. Options:
#'   \itemize{
#'     \item `"none"`: no transformation
#'     \item `"log2"`: log base 2 transformation
#'     \item `"x_1_3"`: cube-root transformation
#'     \item `"sqrt"`: square-root transformation
#'     \item `"boxcox_best"`: data-driven Box-Cox approximation with heuristic labeling
#'   }
#'
#' @return A `MultiAssayExperiment` object with transformed exposures in `colData`,
#' and transformation details stored in:
#' \itemize{
#'   \item `metadata(expomicset)$quality_control$transformation$norm_df`: Shapiro-Wilk test results
#'   \item `metadata(expomicset)$quality_control$transformation$norm_summary`: Summary of normality
#'   \item `metadata(expomicset)$codebook`: Updated with transformation info per variable
#'   \item `metadata(expomicset)$summary$steps`: Updated with step record
#' }
#'
#' @details
#' For `transform_method = "boxcox_best"`, the function automatically shifts values to be
#' strictly positive and chooses from a discrete set of transformations (e.g., `1/x`, `log(x)`, `x^2`)
#' based on estimated Box-Cox lambda. Each variable may receive a different transformation.
#'
#' @examples
#' \dontrun{
#' expomicset <- transform_exposure(expomicset, transform_method = "boxcox_best")
#' }
#'
#' @seealso \code{\link[MASS]{boxcox}}, \code{\link[stats]{shapiro.test}}
#'
#' @export
transform_exposure <- function(
    expomicset,
    exposure_cols = NULL,
    transform_method = "boxcox_best") {

  boxcox_best <- function(var) {
    # Shift to make values strictly positive
    # code credit: https://academic.oup.com/bioinformatics/article/35/24/5344/5523848#supplementary-data
    shift_val <- 0
    if (any(var <= 0, na.rm = TRUE)) {
      shift_val <- abs(min(var, na.rm = TRUE)) + 0.00001
      var <- var + shift_val
    }

    bc <- MASS::boxcox(var ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]

    if (lambda < -1.5) {
      transformed <- 1 / (var^2)
      label <- "1/x^2"
    } else if (lambda < -0.75) {
      transformed <- 1 / var
      label <- "1/x"
    } else if (lambda < -0.25) {
      transformed <- 1 / sqrt(var)
      label <- "1/sqrt(x)"
    } else if (lambda < 0.25) {
      transformed <- log(var)
      label <- "log(x)"
    } else if (lambda < 0.75) {
      transformed <- sqrt(var)
      label <- "sqrt(x)"
    } else if (lambda < 1.5) {
      transformed <- var
      label <- "identity"
    } else {
      transformed <- var^2
      label <- "x^2"
    }

    return(list(values = transformed, label = label))
  }


  apply_transformation <- function(data, columns, method) {
    data <- as.data.frame(data)
    if (method == "none") {
      return(data)
    } else if (method == "log2") {
      return(dplyr::mutate(data, dplyr::across(dplyr::all_of(columns), ~ log2(. + abs(min(.)) + 1))))
    } else if (method == "x_1_3") {
      return(dplyr::mutate(data, dplyr::across(dplyr::all_of(columns), ~ (. + abs(min(.)))^(1/3))))
    } else if (method == "sqrt") {
      return(dplyr::mutate(data, dplyr::across(dplyr::all_of(columns), ~ sqrt(. + abs(min(.))))))
    } else if (method == "boxcox_best") {
      # transformed <- lapply(columns, function(col) boxcox_best(data[[col]]))
      # names(transformed) <- columns
      # return(cbind(data[, setdiff(names(data), columns)], as.data.frame(transformed)))
      # --- Beginning of Test ---
      result <- lapply(columns, function(col) boxcox_best(data[[col]]))
      transformed <- lapply(result, `[[`, "values")
      labels <- sapply(result, `[[`, "label")
      names(transformed) <- columns
      transformed_df <- as.data.frame(transformed)
      attr(transformed_df, "labels") <- labels  # attach to the data.frame AFTER creation
      return(transformed_df)
      # --- End of Test ---
    } else {
      stop("Unsupported transformation method: ", method)
    }
  }

  col_data <- as.data.frame(MultiAssayExperiment::colData(expomicset))
  numeric_cols <- names(col_data)[sapply(col_data, is.numeric)]
  variables_to_transform <- if (!is.null(exposure_cols)) {
    exposure_cols
  } else {
    MultiAssayExperiment::metadata(expomicset)$quality_control$normality$norm_df$exposure
  }
  variables_to_transform <- intersect(variables_to_transform, numeric_cols)

  if (length(variables_to_transform) == 0) {
    stop("No numeric exposure variables to transform.")
  }

  message("Applying the ", transform_method, " transformation.")
  numeric_data <- col_data[, variables_to_transform, drop = FALSE]
  non_numeric_data <- col_data[, setdiff(names(col_data), variables_to_transform)]

  transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)

  if ( transform_method  == "boxcox_best") {
    transformation_info <- data.frame(
      variable = colnames(transformed_numeric),
      transformation_applied = attr(transformed_numeric, "labels"),
      stringsAsFactors = FALSE
    )
  } else{
    transformation_info <- data.frame(
      variable = variables_to_transform,
      transformation_applied = transform_method,
      stringsAsFactors = FALSE
    )
  }

  updated_col_data <- dplyr::bind_cols(non_numeric_data, transformed_numeric)
  MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)

  norm_results <- transformed_numeric |>
    dplyr::select_if(~ length(unique(na.omit(.))) >= 3) |>
    apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
    (\(x) do.call(rbind, x))() |>
    dplyr::mutate(exposure = colnames(transformed_numeric))

  norm_summary <- norm_results |>
    dplyr::summarise(
      "Normal" = sum(p.value > 0.05),
      "Not Normal" = sum(p.value <= 0.05)
    ) |>
    t() |> as.data.frame() |> tibble::rownames_to_column("var") |>
    setNames(c("var", "value"))

  MultiAssayExperiment::metadata(expomicset)$quality_control$transformation <- list(
    norm_df = norm_results,
    norm_summary = norm_summary
  )

  # if (transform_method == "best") {
  #   message("Evaluating the best transformation method including Box-Cox.")
  #
  #   non_numeric_data <- col_data[, setdiff(names(col_data), variables_to_transform)]
  #   numeric_data <- col_data[, variables_to_transform, drop = FALSE]
  #
  #   transformations <- list(
  #     none_trans = apply_transformation(numeric_data, variables_to_transform, "none"),
  #     log_trans = apply_transformation(numeric_data, variables_to_transform, "log2"),
  #     x_1_3_trans = apply_transformation(numeric_data, variables_to_transform, "x_1_3"),
  #     sqrt_trans = apply_transformation(numeric_data, variables_to_transform, "sqrt"),
  #     boxcox_trans = apply_transformation(numeric_data, variables_to_transform, "boxcox_best")
  #   )
  #
  #   norm_results <- lapply(names(transformations), function(name) {
  #     transformed <- transformations[[name]]
  #     transformed |>
  #       apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
  #       (\(x) do.call(rbind, x))() |>
  #       dplyr::mutate(transformation = name,
  #                     exposure = colnames(transformed))
  #   }) |>
  #     dplyr::bind_rows()
  #
  #   norm_summary <- norm_results |>
  #     dplyr::group_by(transformation) |>
  #     dplyr::summarise(
  #       n = dplyr::n(),
  #       normal = sum(p.value > 0.05),
  #       not_normal = sum(p.value <= 0.05),
  #       percent_normal = normal / n
  #     ) |>
  #     dplyr::arrange(dplyr::desc(percent_normal))
  #
  #   best_transformation <- norm_summary$transformation[1]
  #   message("Using the ", best_transformation, " transformation.")
  #
  #   transformed_numeric <- transformations[[best_transformation]]
  #   updated_col_data <- dplyr::bind_cols(non_numeric_data, transformed_numeric)
  #   MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)
  #
  #   MultiAssayExperiment::metadata(expomicset)$quality_control$transformation <- list(
  #     norm_df = norm_results |> dplyr::filter(transformation == best_transformation),
  #     norm_summary = norm_summary
  #   )
  #
  # } else {
  #   message("Applying the ", transform_method, " transformation.")
  #   numeric_data <- col_data[, variables_to_transform, drop = FALSE]
  #   non_numeric_data <- col_data[, setdiff(names(col_data), variables_to_transform)]
  #
  #   transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)
  #
  #   if ( transform_method  == "boxcox_best") {
  #     transformation_info <- data.frame(
  #       variable = colnames(transformed_numeric),
  #       transformation_applied = attr(transformed_numeric, "labels"),
  #       stringsAsFactors = FALSE
  #     )
  #   } else{
  #     transformation_info <- data.frame(
  #       variable = variables_to_transform,
  #       transformation_applied = transform_method,
  #       stringsAsFactors = FALSE
  #     )
  #   }
  #
  #   updated_col_data <- dplyr::bind_cols(non_numeric_data, transformed_numeric)
  #   MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)
  #
  #   norm_results <- transformed_numeric |>
  #     dplyr::select_if(~ length(unique(na.omit(.))) >= 3) |>
  #     apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
  #     (\(x) do.call(rbind, x))() |>
  #     dplyr::mutate(exposure = colnames(transformed_numeric))
  #
  #   norm_summary <- norm_results |>
  #     dplyr::summarise(
  #       "Normal" = sum(p.value > 0.05),
  #       "Not Normal" = sum(p.value <= 0.05)
  #     ) |>
  #     t() |> as.data.frame() |> tibble::rownames_to_column("var") |>
  #     setNames(c("var", "value"))
  #
  #   MultiAssayExperiment::metadata(expomicset)$quality_control$transformation <- list(
  #     norm_df = norm_results,
  #     norm_summary = norm_summary
  #   )
  # }

  # Add step record
  n_transformed <- length(variables_to_transform)

  # Pull normality results
  norm_df <- MultiAssayExperiment::metadata(expomicset)$quality_control$transformation$norm_df

  n_normal <- sum(norm_df$p.value > 0.05, na.rm = TRUE)
  p_normal <- round(n_normal / n_transformed * 100, 1)

  step_record <- list(
    transform_exposure = list(
      timestamp = Sys.time(),
      params = list(
        transform_method = transform_method,
        n_transformed = n_transformed,
        exposures = variables_to_transform
      ),
      notes = paste0(
        "Applied '", transform_method, "' transformation to ", n_transformed, " exposure variables. ",
        n_normal, " passed normality (Shapiro-Wilk p > 0.05, ", p_normal, "%)."
      )
    )
  )

  MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$summary$steps,
    step_record
  )

  # --- Testing this chunk -----

  # Merge with existing codebook or initialize
  existing_codebook <- MultiAssayExperiment::metadata(expomicset)$codebook

  # Create a new data frame for the updated variable information
  updated_codebook <- existing_codebook |>
    full_join(
      transformation_info,
      by = "variable"
    )

  # Update the codebook in the metadata
  MultiAssayExperiment::metadata(expomicset)$codebook <- updated_codebook


  # --- End of Test -------

  return(expomicset)
}


# transform_exposure <- function(
    #     expomicset,
#     exposure_cols=NULL,
#     transform_method = "best") {
#
#   if(!is.null(exposure_cols)){
#     variables_to_transform <- exposure_cols
#   }else{
#     # Extract variables to transform from normality metadata
#     variables_to_transform <- MultiAssayExperiment::metadata(expomicset)$normality$norm_df$exposure
#   }
#
#   if (is.null(variables_to_transform) || length(variables_to_transform) == 0) {
#     stop("Please run `check_normality` before transforming the exposure data.")
#   }
#
#
#   # Helper function to apply transformations
#   apply_transformation <- function(data, columns, method) {
#     if (method == "none") {
#       data
#     } else if (method == "log2") {
#
#       # log2 transformation with minimum value adjustment
#       data |>
#         dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ log2(. + abs(min(.)) + 1)))
#
#     } else if (method == "x_1_3") {
#
#       # x^(1/3) transformation with minimum value adjustment
#       data |>
#         dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ (. + abs(min(.)))^(1/3)))
#
#     } else if (method == "sqrt") {
#
#       # sqrt transformation with minimum value adjustment
#       data |>
#         dplyr::mutate(dplyr::across(dplyr::all_of(columns), ~ sqrt(. + abs(min(.)))))
#
#     } else {
#       stop("Unsupported transformation method: ", method)
#     }
#   }
#
#   # Confirm that variables are numeric
#   col_data <- MultiAssayExperiment::colData(expomicset) |>
#     as.data.frame()
#   numeric_cols <- col_data |>
#     dplyr::select_if(~is.numeric(.)) |>
#     colnames()
#   variables_to_transform <- variables_to_transform[variables_to_transform %in% numeric_cols]
#
#   # Perform transformations based on the method
#   if (transform_method == "best") {
#     message("Evaluating the best transformation method including no transformation.")
#
#     # Separate numeric and non-numeric columns
#     col_data <- MultiAssayExperiment::colData(expomicset) |>
#       as.data.frame()
#     numeric_data <- col_data |>
#       dplyr::select(dplyr::all_of(variables_to_transform))
#     non_numeric_data <- col_data |>
#       dplyr::select(-dplyr::all_of(variables_to_transform))
#
#     # Apply transformations to numeric columns
#     transformations <- list(
#       none_trans = apply_transformation(numeric_data, variables_to_transform, "none"),
#       log_trans = apply_transformation(numeric_data, variables_to_transform, "log2"),
#       x_1_3_trans = apply_transformation(numeric_data, variables_to_transform, "x_1_3"),
#       sqrt_trans = apply_transformation(numeric_data, variables_to_transform, "sqrt")
#     )
#
#     # Evaluate normality for each transformation
#     norm_results <- lapply(names(transformations), function(name) {
#       transformed <- transformations[[name]]
#       transformed |>
#         apply(2, function(x) shapiro.test(x) |>
#                 broom::tidy()) |>
#         (\(x) do.call(rbind, x))() |>
#         dplyr::mutate(transformation = name,
#                exposure = colnames(transformed))
#     }) |>
#       bind_rows()
#
#     # Summarize normality results
#     norm_summary <- norm_results |>
#       dplyr::group_by(transformation) |>
#       dplyr::summarise(
#         n = dplyr::n(),
#         normal = sum(p.value > 0.05),
#         not_normal = sum(p.value <= 0.05),
#         percent_normal = normal / n
#       ) |>
#       dplyr::arrange(dplyr::desc(percent_normal))
#
#     # Select the best transformation
#     best_transformation <- norm_summary$transformation[1]
#     if (best_transformation == "none_trans") {
#       message("No transformation applied.")
#     } else {
#       message("Using the ", best_transformation, " transformation.")
#     }
#
#     # Update colData with the best transformation
#     transformed_numeric <- transformations[[best_transformation]]
#     updated_col_data <- non_numeric_data |>
#       dplyr::bind_cols(transformed_numeric)
#     MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)
#
#     # only save the results for the chosen transformation
#     norm_results <- norm_results |>
#       dplyr::filter(transformation == best_transformation)
#
#     # Save results
#     MultiAssayExperiment::metadata(expomicset)$transformation <- list(
#       norm_df = norm_results,
#       norm_summary = norm_summary
#     )
#
#   } else {
#     # Apply a specific transformation
#     message("Applying the ", transform_method, " transformation.")
#
#     # Separate numeric and non-numeric columns
#     col_data <- MultiAssayExperiment::colData(expomicset) |>
#       as.data.frame()
#     numeric_data <- col_data |>
#       dplyr::select(dplyr::all_of(variables_to_transform))
#     non_numeric_data <- col_data |>
#       dplyr::select(-dplyr::all_of(variables_to_transform))
#
#     # Transform the numeric data
#     transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)
#
#     # Combine with non-numeric columns
#     updated_col_data <- non_numeric_data |>
#       bind_cols(transformed_numeric)
#     MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)
#
#     # Perform Shapiro-Wilk test
#     norm_results <- transformed_numeric |>
#       dplyr::select_if(~ length(unique(na.omit(.))) >= 3) |>
#       apply(2, function(x) shapiro.test(x) |>
#               broom::tidy()) |>
#       (\(x) do.call(rbind, x))() |>
#       dplyr::mutate(exposure = colnames(transformed_numeric))
#
#     # Summarize normality results
#     norm_summary <- norm_results |>
#       dplyr::summarise(
#         "Normal" = sum(p.value > 0.05),
#         "Not Normal" = sum(p.value <= 0.05)
#       ) |>
#       t() |>
#       as.data.frame() |>
#       rownames_to_column("var") |>
#       setNames(c("var","value"))
#
#
#
#     # Save results
#     MultiAssayExperiment::metadata(expomicset)$transformation <- list(
#       norm_df = norm_results,
#       norm_summary = norm_summary
#     )
#   }
#
#   return(expomicset)
# }
