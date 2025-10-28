#' Perform Principal Component Analysis (PCA)
#'
#' Runs PCA on the feature and sample spaces of a `MultiAssayExperiment` object,
#' identifying outliers based on Mahalanobis distance.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing omics
#' and exposure data.
#' @param log_trans_exp A boolean value specifying whether to log2
#' transform the exposure data
#' @param log_trans_omics a boolean value specifying whether to log2
#' transform the omics data
#' @param action A character string specifying whether to store
#' (`"add"`) or return (`"get"`) the results. Default is `"add"`.
#'
#' @details
#' This function:
#' - Identifies **common samples** across all assays and exposure data.
#' - Performs **PCA on features** (transformed and standardized).
#' - Performs **PCA on samples** and computes Mahalanobis distance to
#' detect outliers.
#' - **Output Handling**:
#'   - `"add"`: Stores results in `metadata(exposomicset)$pca` and
#'   updates `colData` with PCs.
#'   - `"get"`: Returns a list containing the PCA results.
#'
#' @return A `MultiAssayExperiment` object with PCA results added to
#' metadata (if `action = "add"`) or a list with:
#' \item{pca_df}{A tibble of the transformed input data.}
#' \item{pca_feature}{A `prcomp` object containing PCA results for features.}
#' \item{pca_sample}{A `prcomp` object containing PCA results for samples.}
#' \item{outliers}{A character vector of detected sample outliers.}
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run pca
#' mae <- mae |>
#'     run_pca()
#'
#' @export
run_pca <- function(
  exposomicset,
  log_trans_exp = FALSE,
  log_trans_omics = TRUE,
  action = "add"
) {
    # Identify common samples across all data
    message("Identifying common samples.")

    common_samples <- rownames(MultiAssayExperiment::colData(exposomicset))
    for (omics_name in names(MultiAssayExperiment::experiments(exposomicset))) {
        common_samples <- intersect(
            common_samples,
            colnames(MultiAssayExperiment::experiments(
                exposomicset
            )[[omics_name]])
        )
    }

    if (length(common_samples) == 0) {
        stop("No common samples found across exposure and omics data.")
    }

    # Subset colData to common samples
    message("Subsetting exposure data.")

    exposure_data <- MultiAssayExperiment::colData(
        exposomicset
    )[common_samples, ] |>
        as.data.frame() |>
        dplyr::select_if(is.numeric) |>
        t() |>
        as.data.frame()
    exposure_data <- transform(exposure_data, category = "exposure")

    # Remove exposures with any NA values
    exposure_data <- exposure_data[rowSums(is.na(exposure_data)) == 0, , drop = FALSE]

    if (log_trans_exp) {
        exposure_data <- exposure_data |>
            dplyr::mutate(
                dplyr::across(dplyr::where(is.numeric), ~ .safe_log2(. + abs(min(.)) + 1))
            )
    }

    # Subset omics data to common samples
    message("Subsetting omics data.")

    omics_data <- lapply(
        names(MultiAssayExperiment::experiments(exposomicset)),
        function(omics_name) {
            SummarizedExperiment::assays(
                MultiAssayExperiment::experiments(exposomicset)[[omics_name]]
            )[[1]][
                , common_samples,
                drop = FALSE
            ] |>
                as.data.frame() |>
                transform(category = omics_name)
        }
    ) |>
        dplyr::bind_rows()

    if (log_trans_omics) {
        omics_data <- omics_data |>
            dplyr::mutate(
                dplyr::across(dplyr::where(is.numeric), ~ .safe_log2(. + abs(min(.)) + 1))
            )
    }

    # Combine datasets
    dat <- exposure_data |>
        dplyr::bind_rows(omics_data)
    dat <- transform(dat, id = rownames(dat))

    # Remove zero-variance columns
    feature_data <- dat |>
        dplyr::select(-c(category, id)) |>
        dplyr::select(where(\(x) var(x, na.rm = TRUE) > 0))

    # PCA analysis: feature space
    message("Performing PCA on Feature Space.")
    pca_feature <- prcomp(feature_data, center = TRUE, scale. = TRUE)


    # PCA analysis: sample space
    message("Performing PCA on Sample Space.")
    sample_data <- dat |>
        dplyr::select(-c(category, id)) |>
        t() |>
        as.data.frame()

    sample_vars <- apply(sample_data, 2, \(x) var(x, na.rm = TRUE))
    sample_data <- sample_data |>
        dplyr::select(names(sample_vars[sample_vars > 0]))

    pca_sample <- prcomp(sample_data, center = TRUE, scale. = TRUE)

    # Compute Mahalanobis distance
    pca_coords <- pca_sample$x[, c("PC1", "PC2")] # Adjust if using more PCs
    cov_mat <- cov(pca_coords) # Covariance matrix
    inv_cov <- solve(cov_mat) # Inverse covariance
    center <- colMeans(pca_coords) # Mean vector

    mahal_dist <- apply(pca_coords, 1, function(x) {
        sqrt(t(x - center) %*% inv_cov %*% (x - center))
    })

    # Determine threshold (Chi-square with 2 degrees of freedom, p < 0.01)
    threshold <- qchisq(0.99, df = 2)
    outliers <- which(mahal_dist^2 > threshold)

    # Print outlier sample names
    if (length(outliers) > 0) {
        message(paste(
            "Outliers detected:",
            paste(rownames(pca_sample$x)[outliers],
                collapse = ", "
            )
        ))
    } else {
        message("No outliers detected.")
    }

    # add in PCs for samples to the colData
    col_data <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        (\(df){
            df$id_to_map <- rownames(df)
            df
        })() |>
        dplyr::left_join(
            pca_sample$x |>
                as.data.frame() |>
                (\(df){
                    df$id_to_map <- rownames(df)
                    df
                })(),
            by = "id_to_map"
        ) |>
        tibble::column_to_rownames("id_to_map")

    MultiAssayExperiment::colData(exposomicset) <- col_data

    if (action == "add") {
        # Store results
        MultiAssayExperiment::metadata(exposomicset)$quality_control$pca <- list(
            pca_df = tibble(dat),
            pca_feature = pca_feature,
            pca_sample = pca_sample,
            outliers = rownames(pca_sample$x)[outliers]
        )

        # Add analysis steps taken to metadata
        step_record <- list(
            run_pca = list(
                timestamp = Sys.time(),
                params = list(),
                notes = paste("Outliers: ",
                    paste(rownames(pca_sample$x)[outliers], collapse = ", "),
                    collapse = ""
                )
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else if (action == "get") {
        return(list(
            pca_df = tibble(dat),
            pca_feature = pca_feature,
            pca_sample = pca_sample,
            outliers = rownames(pca_sample$x)[outliers]
        ))
    } else {
        stop("Invalid action. Use 'add' or 'get'.")
    }
}


#' Safe Log2 Transformation
#'
#' This helper function ensures that negative, infinite, or `NaN` values
#' do not cause errors during log transformation. It converts all
#' negative values to zero, replaces `Inf`/`NaN` with `NA`, and performs
#' `log2(x + 1)` on the adjusted input.
#'
#' @param x A numeric vector to be log2-transformed.
#' @return A numeric vector of the same length as `x`, with safely
#' log2-transformed values.
#'
#' @noRd
#' @keywords internal
.safe_log2 <- function(x) {
    x <- as.numeric(x)
    x <- ifelse(is.infinite(x) | is.nan(x), NA, x)
    log2(pmax(x, 0) + 1)
}
