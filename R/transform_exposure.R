#' Transform Exposure Variables in a MultiAssayExperiment
#'
#' Applies various transformation methods to numeric exposure variables in the
#' `colData` of a `MultiAssayExperiment` object to improve normality.
#'
#' Supported transformations include:
#' \itemize{
#'   \item \code{"none"}: No transformation
#'   \item \code{"log2"}: Log base 2 with minimum value shift
#'   \item \code{"x_1_3"}: Cube root with minimum value shift
#'   \item \code{"sqrt"}: Square root with minimum value shift
#'   \item \code{"boxcox_best"}: Box-Cox power transformation based on optimal lambda
#'   \item \code{"best"}: Automatically selects the transformation method that best
#'         normalizes the data (based on Shapiro-Wilk test)
#' }
#'
#' The \code{"boxcox_best"} transformation is based on the Box-Cox method from the
#' \code{MASS} package. Variables are shifted to ensure strictly positive values before
#' transformation. Lambda intervals are interpreted as follows:
#' \itemize{
#'   \item \eqn{\lambda < -1.5}: \eqn{1 / x^2}
#'   \item \eqn{-1.5 \leq \lambda < -0.75}: \eqn{1 / x}
#'   \item \eqn{-0.75 \leq \lambda < -0.25}: \eqn{1 / sqrt(x)}
#'   \item \eqn{-0.25 \leq \lambda < 0.25}: \eqn{log(x)}
#'   \item \eqn{0.25 \leq \lambda < 0.75}: \eqn{sqrt(x)}
#'   \item \eqn{0.75 \leq \lambda < 1.5}: identity
#'   \item \eqn{\lambda \geq 1.5}: \eqn{x^2}
#' }
#'
#' @param expomicset A \code{MultiAssayExperiment} object.
#' @param cols_of_interest A character vector of column names from \code{colData(expomicset)} to transform.
#'   If \code{NULL}, the function will attempt to use metadata from \code{check_normality()}.
#' @param transform_method One of \code{"none"}, \code{"log2"}, \code{"x_1_3"}, \code{"sqrt"},
#'   \code{"boxcox_best"}, or \code{"best"} (default). \code{"best"} evaluates all available
#'   methods and chooses the one that most improves normality.
#'
#' @return A \code{MultiAssayExperiment} object with updated \code{colData} and transformation
#'   results saved in \code{metadata(expomicset)$transformation}.
#'
#' @importFrom MultiAssayExperiment colData metadata
#' @importFrom dplyr mutate across all_of bind_cols select select_if group_by summarise n arrange desc
#' @importFrom broom tidy
#' @importFrom MASS boxcox
#' @importFrom tibble rownames_to_column
#' @export
#'
#' @examples
#' \dontrun{
#' transformed_mae <- transform_exposure(my_mae,
#'                                       cols_of_interest = c("pm25", "no2"),
#'                                       transform_method = "best")
#' }

transform_exposure <- function(
    expomicset,
    cols_of_interest = NULL,
    transform_method = "best") {

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
    } else if (lambda < -0.75) {
      transformed <- 1 / var
    } else if (lambda < -0.25) {
      transformed <- 1 / sqrt(var)
    } else if (lambda < 0.25) {
      transformed <- log(var)
    } else if (lambda < 0.75) {
      transformed <- sqrt(var)
    } else if (lambda < 1.5) {
      transformed <- var
    } else {
      transformed <- var^2
    }

    return(transformed)
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
      transformed <- lapply(columns, function(col) boxcox_best(data[[col]]))
      names(transformed) <- columns
      return(cbind(data[, setdiff(names(data), columns)], as.data.frame(transformed)))
    } else {
      stop("Unsupported transformation method: ", method)
    }
  }

  col_data <- as.data.frame(MultiAssayExperiment::colData(expomicset))
  numeric_cols <- names(col_data)[sapply(col_data, is.numeric)]
  variables_to_transform <- if (!is.null(cols_of_interest)) {
    cols_of_interest
  } else {
    MultiAssayExperiment::metadata(expomicset)$normality$norm_df$exposure
  }
  variables_to_transform <- intersect(variables_to_transform, numeric_cols)

  if (length(variables_to_transform) == 0) {
    stop("No numeric exposure variables to transform.")
  }

  if (transform_method == "best") {
    message("Evaluating the best transformation method including Box-Cox.")

    non_numeric_data <- col_data[, setdiff(names(col_data), variables_to_transform)]
    numeric_data <- col_data[, variables_to_transform, drop = FALSE]

    transformations <- list(
      none_trans = apply_transformation(numeric_data, variables_to_transform, "none"),
      log_trans = apply_transformation(numeric_data, variables_to_transform, "log2"),
      x_1_3_trans = apply_transformation(numeric_data, variables_to_transform, "x_1_3"),
      sqrt_trans = apply_transformation(numeric_data, variables_to_transform, "sqrt"),
      boxcox_trans = apply_transformation(numeric_data, variables_to_transform, "boxcox_best")
    )

    norm_results <- lapply(names(transformations), function(name) {
      transformed <- transformations[[name]]
      transformed |>
        apply(2, function(x) shapiro.test(x) |> broom::tidy()) |>
        (\(x) do.call(rbind, x))() |>
        dplyr::mutate(transformation = name,
                      exposure = colnames(transformed))
    }) |>
      dplyr::bind_rows()

    norm_summary <- norm_results |>
      dplyr::group_by(transformation) |>
      dplyr::summarise(
        n = dplyr::n(),
        normal = sum(p.value > 0.05),
        not_normal = sum(p.value <= 0.05),
        percent_normal = normal / n
      ) |>
      dplyr::arrange(dplyr::desc(percent_normal))

    best_transformation <- norm_summary$transformation[1]
    message("Using the ", best_transformation, " transformation.")

    transformed_numeric <- transformations[[best_transformation]]
    updated_col_data <- dplyr::bind_cols(non_numeric_data, transformed_numeric)
    MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(updated_col_data)

    MultiAssayExperiment::metadata(expomicset)$transformation <- list(
      norm_df = norm_results |> dplyr::filter(transformation == best_transformation),
      norm_summary = norm_summary
    )

  } else {
    message("Applying the ", transform_method, " transformation.")
    numeric_data <- col_data[, variables_to_transform, drop = FALSE]
    non_numeric_data <- col_data[, setdiff(names(col_data), variables_to_transform)]

    transformed_numeric <- apply_transformation(numeric_data, variables_to_transform, transform_method)
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

    MultiAssayExperiment::metadata(expomicset)$transformation <- list(
      norm_df = norm_results,
      norm_summary = norm_summary
    )
  }

  # Add analysis steps taken to metadata
  MultiAssayExperiment::metadata(expomicset)$steps <- c(
    MultiAssayExperiment::metadata(expomicset)$steps,
    "transform_exposure"
  )

  return(expomicset)
}


# transform_exposure <- function(
#     expomicset,
#     cols_of_interest=NULL,
#     transform_method = "best") {
#
#   if(!is.null(cols_of_interest)){
#     variables_to_transform <- cols_of_interest
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
