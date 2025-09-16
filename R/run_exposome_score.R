#' Compute Composite Exposome Scores
#'
#' Calculates a summary exposome score per sample using one of several methods
#' including mean, sum, median, PCA (first principal component),
#' IRT (Item Response Theory), quantile binning, or row-wise variance.
#' The resulting score is added to the `colData` of the
#' `MultiAssayExperiment` object.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing exposure
#' data in its `colData`.
#' @param score_type Character. The method used to compute the score.
#' Options are:
#'   `"mean"`, `"sum"`, `"median"`, `"pca"`, `"irt"`, `"quantile"`, `"var"`.
#' @param exposure_cols Optional character vector. Specific exposure
#'  column names to include. If `NULL`, all numeric columns are used.
#' @param scale Logical. Whether to scale the exposures before computing
#' the score. Default is `TRUE`.
#' @param score_column_name Optional name for the resulting score column.
#' If `NULL`, an automatic name is used (e.g., `"exposome_score_pca"`).
#'
#' @return A `MultiAssayExperiment` object with the exposome score added
#'  to `colData()`.
#'
#' @details
#' - `"pca"` uses the first principal component from `prcomp()`.
#' - `"irt"` uses the `mirt` package to fit a graded response model to
#' discretized exposures.
#' - `"quantile"` assigns decile bins (1-10) to each variable and sums
#' them row-wise.
#' - `"var"` computes the row-wise variance across exposures.
#'
#' @examples
#' # create the example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # create the air pollution score
#' mae <- run_exposome_score(
#'     mae,
#'     score_type = "pca",
#'     exposure_cols = c("exposure_pm25", "exposure_no2"),
#'     scale = TRUE,
#'     score_column_name = "air_pollution_score"
#' )
#'
#' @export
run_exposome_score <- function(
    exposomicset,
    score_type,
    exposure_cols = NULL,
    scale = TRUE,
    score_column_name = NULL) {
    # Extract and preprocess colData
    message("Extracting exposure data...")
    data <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        dplyr::select_if(is.numeric)

    # Unselect PCs
    data <- data |>
        dplyr::select(-dplyr::starts_with("PC"))

    # Check if exposure_cols is provided
    if (!is.null(exposure_cols)) {
        exposure_cols <- intersect(exposure_cols, colnames(data))

        data <- data |>
            dplyr::select(all_of(exposure_cols))
    }

    # Scale the data if user says so
    if (scale) {
        data <- data |>
            dplyr::mutate_all(~ scale(.))
    }

    if (score_type == "mean") {
        # Calculate mean exposome score
        message("Calculating mean exposure scores...")

        scores <- data |>
            dplyr::mutate(exposome_score_mean = rowMeans(
                dplyr::across(dplyr::everything()),
                na.rm = TRUE
            )) |>
            dplyr::select(exposome_score_mean)
    } else if (score_type == "sum") {
        # Calculate sum exposome score
        message("Calculating sum exposure scores...")

        scores <- data |>
            dplyr::mutate(exposome_score_sum = rowSums(
                dplyr::across(dplyr::everything()),
                na.rm = TRUE
            )) |>
            dplyr::select(exposome_score_sum)
    } else if (score_type == "median") {
        # Calculate median exposome score
        message("Calculating median exposure scores...")
        .check_suggested("matrixStats")
        scores <- data |>
            dplyr::mutate(exposome_score_median = matrixStats::rowMedians(
                as.matrix(
                    dplyr::pick(
                        dplyr::everything()
                    )
                ),
                na.rm = TRUE
            )) |>
            dplyr::select(exposome_score_median)
    } else if (score_type == "pca") {
        # Calculate PCA exposome score
        message("Calculating PCA exposure scores...")

        pca_result <- prcomp(data, center = TRUE)
        scores <- data.frame(pca_result$x[, 1])
        colnames(scores) <- c("exposome_score_pca")
    } else if (score_type == "irt") {
        # Calculate IRT exposome score
        message("Calculating IRT exposure scores...")
        .check_suggested("mirt")

        exposures_ordinal <- apply(data, 2, function(x) dplyr::ntile(x, 10))
        model <- mirt::mirt(exposures_ordinal, 1, itemtype = "graded")
        scores <- mirt::fscores(model, full.scores.SE = FALSE)
        scores <- scores |>
            as.data.frame() |>
            setNames("exposome_score_irt")
        rownames(scores) <- rownames(data)
    } else if (score_type == "quantile") {
        # Calculate quantile exposome score
        message("Calculating quantile exposure scores...")

        quantile_data <- data |>
            dplyr::mutate(
                dplyr::across(
                    dplyr::everything(),
                    ~ dplyr::ntile(., 10)
                )
            )

        scores <- quantile_data |>
            dplyr::mutate(exposome_score_quantile = rowSums(
                dplyr::across(dplyr::everything()),
                na.rm = TRUE
            )) |>
            dplyr::select(exposome_score_quantile)
    } else if (score_type == "var") {
        # Calculate variance exposome score
        message("Calculating variance exposure scores...")
        .check_suggested("matrixStats")
        scores <- data |>
            dplyr::mutate(exposome_score_var = matrixStats::rowVars(
                as.matrix(
                    dplyr::pick(
                        dplyr::everything()
                    )
                ),
                na.rm = TRUE
            )) |>
            dplyr::select(exposome_score_var)
    } else {
        stop("Choose either 'sum', mean', 'median', 'pca', 'irt', or 'quantile'.")
    }

    if (!is.null(score_column_name)) {
        colnames(scores) <- score_column_name
    }

    # Add scores to the MultiAssayExperiment object
    updated_col_data <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        # if name exists in df - overwrite the score column name
        dplyr::select(-dplyr::any_of(colnames(scores))) |>
        tibble::rownames_to_column("id") |>
        dplyr::left_join(
            scores |>
                tibble::rownames_to_column("id"),
            by = "id"
        ) |>
        tibble::column_to_rownames("id")

    # Update the colData of the MultiAssayExperiment object
    MultiAssayExperiment::colData(exposomicset) <- S4Vectors::DataFrame(
        updated_col_data
    )

    # Add summary step record
    # Create step content
    step_content <- list(
        timestamp = Sys.time(),
        params = list(
            score_type = score_type,
            scaled = scale,
            exposure_cols = exposure_cols,
            score_column_name = score_column_name
        ),
        notes = paste0(
            "Exposome score computed using method: '",
            score_type,
            "'"
        )
    )

    # Assign dynamic name to the step
    step_record <- setNames(
        list(step_content),
        paste0(
            "run_exposome_score_",
            score_column_name
        )
    )

    MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
        MultiAssayExperiment::metadata(exposomicset)$summary$steps,
        step_record
    )

    return(exposomicset)
}
