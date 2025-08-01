#' Transform Exposure Variables for Normality
#'
#' Applies a transformation to selected numeric exposure variables in the
#' `colData` of a `MultiAssayExperiment` to improve their normality
#' (e.g., log, Box-Cox, sqrt). Transformation results and normality
#' statistics are stored in metadata for tracking.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure
#' variables in `colData`.
#' @param exposure_cols Optional character vector of exposure variable
#' names to transform.
#'   If `NULL`, uses exposures found in `metadata(expomicset)$quality_control$normality$norm_df$exposure`.
#' @param transform_method Character. Transformation method to apply. Options:
#'   \itemize{
#'     \item `"none"`: no transformation
#'     \item `"log2"`: log base 2 transformation
#'     \item `"x_1_3"`: cube-root transformation
#'     \item `"sqrt"`: square-root transformation
#'     \item `"boxcox_best"`: data-driven Box-Cox approximation with
#'     heuristic labeling
#'   }
#'
#' @return A `MultiAssayExperiment` object with transformed exposures
#' in `colData`, and transformation details stored in:
#' \itemize{
#'   \item `metadata(expomicset)$quality_control$transformation$norm_df`:
#'    Shapiro-Wilk test results
#'   \item `metadata(expomicset)$quality_control$transformation$norm_summary`:
#'   Summary of normality
#'   \item `metadata(expomicset)$codebook`: Updated with transformation info
#'    per variable
#'   \item `metadata(expomicset)$summary$steps`: Updated with step record
#' }
#'
#' @details
#' For `transform_method = "boxcox_best"`, the function automatically shifts
#' values to be strictly positive and chooses from a discrete set of
#' transformations (e.g., `1/x`, `log(x)`, `x^2`)
#' based on estimated Box-Cox lambda. Each variable may receive a
#' different transformation.
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # Test for normality
#' mae <- mae |>
#'     run_normality_check() |>
#'     transform_exposure(
#'         exposure_cols = c("age", "bmi", "exposure_pm25"),
#'         transform_method = "boxcox_best"
#'     )
#'
#' @seealso \code{\link[MASS]{boxcox}}, \code{\link[stats]{shapiro.test}}
#'
#' @export
transform_exposure <- function(
    expomicset,
    exposure_cols = NULL,
    transform_method = "boxcox_best") {
    col_data <- as.data.frame(MultiAssayExperiment::colData(expomicset))
    # numeric_cols <- names(col_data)[sapply(col_data, is.numeric)]
    numeric_cols <- names(col_data)[
        vapply(col_data, is.numeric, FUN.VALUE = logical(1))
    ]
    variables_to_transform <- if (!is.null(exposure_cols)) {
        exposure_cols
    } else {
        MultiAssayExperiment::metadata(
            expomicset
        )$quality_control$normality$norm_df$exposure
    }
    variables_to_transform <- intersect(variables_to_transform, numeric_cols)

    if (length(variables_to_transform) == 0) {
        stop("No numeric exposure variables to transform.")
    }

    message("Applying the ", transform_method, " transformation.")
    numeric_data <- col_data[, variables_to_transform, drop = FALSE]
    non_numeric_data <- col_data[
        , setdiff(names(col_data), variables_to_transform)
    ]

    transformed_numeric <- .apply_transformation(
        numeric_data,
        variables_to_transform,
        transform_method
    )

    if (transform_method == "boxcox_best") {
        transformation_info <- data.frame(
            variable = colnames(transformed_numeric),
            transformation_applied = attr(transformed_numeric, "labels"),
            stringsAsFactors = FALSE
        )
    } else {
        transformation_info <- data.frame(
            variable = variables_to_transform,
            transformation_applied = transform_method,
            stringsAsFactors = FALSE
        )
    }

    updated_col_data <- dplyr::bind_cols(
        non_numeric_data,
        transformed_numeric
    )
    MultiAssayExperiment::colData(expomicset) <- S4Vectors::DataFrame(
        updated_col_data
    )

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
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column("var") |>
        setNames(c("var", "value"))

    all_metadata <- MultiAssayExperiment::metadata(expomicset)
    all_metadata$quality_control$transformation <- list(
        norm_df = norm_results,
        norm_summary = norm_summary
    )
    MultiAssayExperiment::metadata(expomicset) <- all_metadata

    # Add step record
    n_transformed <- length(variables_to_transform)

    # Pull normality results
    all_metadata <- MultiAssayExperiment::metadata(expomicset)
    norm_df <- all_metadata$quality_control$transformation$norm_df

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
                "Applied '",
                transform_method,
                "' transformation to ",
                n_transformed,
                " exposure variables. ",
                n_normal,
                " passed normality (Shapiro-Wilk p > 0.05, ",
                p_normal,
                "%)."
            )
        )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
        MultiAssayExperiment::metadata(expomicset)$summary$steps,
        step_record
    )

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

    return(expomicset)
}

# --- Box-Cox Best Function ------
#' Box-Cox based transformation selector
#'
#' This internal function applies a Box-Cox analysis to determine the most
#' appropriate transformation for a numeric vector. Depending on the estimated
#' Box-Cox lambda, it applies a power, inverse, log, or identity transform.
#'
#' @param var A numeric vector to transform.
#'
#' @return A list with two elements:
#'   - `values`: the transformed vector
#'   - `label`: the transformation applied as a string
#' @importFrom MASS boxcox
#' @keywords internal
#' @noRd
.boxcox_best <- function(var) {
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

# --- Apply Transformation Function -----
#' Apply a transformation to selected columns
#'
#' This internal function applies a variety of transformations (log, sqrt,
#' Box-Cox, etc.) to specific columns in a data frame.
#'
#' @param data A data frame or matrix containing numeric columns.
#' @param columns Character vector of column names to transform.
#' @param method Transformation method. Supported options:
#'   - `"none"`: no transformation
#'   - `"log2"`: log2(x + |min(x)| + 1)
#'   - `"x_1_3"`: cube root transform (x + |min(x)|)^(1/3)
#'   - `"sqrt"`: square root transform (x + |min(x)|)
#'   - `"boxcox_best"`: adaptive Box-Cox transformation using
#'   \code{.boxcox_best()}
#'
#' @return A data frame with transformed columns. If `method = "boxcox_best"`,
#'   an attribute `"labels"` is attached containing the transformation labels
#'   used for each column.
#'
#' @keywords internal
#' @importFrom dplyr mutate across all_of
#' @noRd
.apply_transformation <- function(data, columns, method) {
    data <- as.data.frame(data)
    if (method == "none") {
        return(data)
    } else if (method == "log2") {
        return(dplyr::mutate(
            data,
            dplyr::across(
                dplyr::all_of(columns),
                ~ log2(. + abs(min(.)) + 1)
            )
        ))
    } else if (method == "x_1_3") {
        return(dplyr::mutate(
            data,
            dplyr::across(
                dplyr::all_of(columns),
                ~ (. + abs(min(.)))^(1 / 3)
            )
        ))
    } else if (method == "sqrt") {
        return(dplyr::mutate(
            data,
            dplyr::across(
                dplyr::all_of(columns),
                ~ sqrt(. + abs(min(.)))
            )
        ))
    } else if (method == "boxcox_best") {
        result <- lapply(columns, function(col) .boxcox_best(data[[col]]))
        transformed <- lapply(result, `[[`, "values")

        labels <- vapply(result, `[[`, "label", FUN.VALUE = character(1))

        names(transformed) <- columns
        transformed_df <- as.data.frame(transformed)

        # attach labels to the data.frame AFTER creation
        attr(transformed_df, "labels") <- labels
        return(transformed_df)
    } else {
        stop("Unsupported transformation method: ", method)
    }
}
