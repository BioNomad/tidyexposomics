#' Perform Principal Component Analysis (PCA)
#'
#' Runs PCA on the feature and sample spaces of a `MultiAssayExperiment` object,
#' identifying outliers based on Mahalanobis distance.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param action A character string specifying whether to store (`"add"`) or return (`"get"`) the results. Default is `"add"`.
#'
#' @details
#' This function:
#' - Identifies **common samples** across all assays and exposure data.
#' - Performs **PCA on features** (transformed and standardized).
#' - Performs **PCA on samples** and computes Mahalanobis distance to detect outliers.
#' - **Output Handling**:
#'   - `"add"`: Stores results in `metadata(expomicset)$pca` and updates `colData` with PCs.
#'   - `"get"`: Returns a list containing the PCA results.
#'
#' @return A `MultiAssayExperiment` object with PCA results added to metadata (if `action = "add"`) or a list with:
#' \item{pca_df}{A tibble of the transformed input data.}
#' \item{pca_feature}{A `prcomp` object containing PCA results for features.}
#' \item{pca_sample}{A `prcomp` object containing PCA results for samples.}
#' \item{outliers}{A character vector of detected sample outliers.}
#'
#' @examples
#' \dontrun{
#' expom <- run_pca(expomicset = expom, action = "add")
#' pca_results <- run_pca(expomicset = expom, action = "get")
#' }
#'
#' @export
run_pca <- function(
    expomicset,
    action="add") {

  # Identify common samples across all data
  message("Identifying common samples...")

  common_samples <- rownames(MultiAssayExperiment::colData(expomicset))
  for (omics_name in names(MultiAssayExperiment::experiments(expomicset))) {
    common_samples <- intersect(common_samples, colnames(MultiAssayExperiment::experiments(expomicset)[[omics_name]]))
  }

  if (length(common_samples) == 0) {
    stop("No common samples found across exposure and omics data.")
  }

  # Subset colData to common samples
  message("Subsetting exposure data...")

  exposure_data <- MultiAssayExperiment::colData(expomicset)[common_samples, ] |>
    as.data.frame() |>
    dplyr::select(where(is.numeric)) |>
    t() |>
    as.data.frame()
  exposure_data <- transform(exposure_data, category = "exposure")

  # Subset omics data to common samples
  message("Subsetting omics data...")

  omics_data <- lapply(
    names(MultiAssayExperiment::experiments(expomicset)),
    function(omics_name) {
      SummarizedExperiment::assays(
        MultiAssayExperiment::experiments(expomicset)[[omics_name]])[[1]][, common_samples, drop = FALSE] |>
        as.data.frame() |>
        transform(category = omics_name)
    }) |>
    dplyr::bind_rows()

  # Combine datasets
  dat <- exposure_data |>
    dplyr::bind_rows(omics_data)
  dat <- transform(dat, id = rownames(dat))

  # Remove zero-variance columns
  feature_data <- dat |>
    dplyr::select(-c(category, id)) |>
    dplyr::select(where(\(x) var(x, na.rm = TRUE) > 0)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ log2(.+abs(min(.))+1)))

  # PCA analysis: feature space
  message("Performing PCA on Feature Space...")
  pca_feature <- prcomp(feature_data, center = TRUE, scale. = TRUE)


  # PCA analysis: sample space
  message("Performing PCA on Sample Space...")
  sample_data <- dat |>
    dplyr::select(-c(category, id)) |>
    t() |>
    as.data.frame()

  sample_vars <- apply(sample_data, 2, \(x) var(x, na.rm = TRUE))
  sample_data <- sample_data |>
    dplyr::select(names(sample_vars[sample_vars > 0])) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ log2(.+abs(min(.))+1)))

  pca_sample <- prcomp(sample_data, center = TRUE, scale. = TRUE)

  # Compute Mahalanobis distance
  pca_coords <- pca_sample$x[, c("PC1", "PC2")]  # Adjust if using more PCs
  cov_mat <- cov(pca_coords)  # Covariance matrix
  inv_cov <- solve(cov_mat)  # Inverse covariance
  center <- colMeans(pca_coords)  # Mean vector

  mahal_dist <- apply(pca_coords, 1, function(x) {
    sqrt(t(x - center) %*% inv_cov %*% (x - center))
  })

  # Determine threshold (Chi-square with 2 degrees of freedom, p < 0.01)
  threshold <- qchisq(0.99, df = 2)
  outliers <- which(mahal_dist^2 > threshold)

  # Print outlier sample names
  if (length(outliers) > 0) {
    print(paste("Outliers detected:",
                paste(rownames(pca_sample$x)[outliers],
                      collapse = ", ")))
  } else {
    print("No outliers detected.")
  }

  # add in PCs for samples to the colData
  col_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    (\(df){df$id_to_map=rownames(df);df})() |>
    dplyr::left_join(pca_sample$x |>
                as.data.frame() |>
                (\(df){df$id_to_map=rownames(df);df})(),
              by="id_to_map") |>
    tibble::column_to_rownames("id_to_map")

  MultiAssayExperiment::colData(expomicset) <- col_data

  if(action=="add"){
    # Store results
    MultiAssayExperiment::metadata(expomicset)$pca <- list(
      pca_df = tibble(dat),
      pca_feature = pca_feature,
      pca_sample = pca_sample,
      outliers = rownames(pca_sample$x)[outliers]
    )

    # Add analysis steps taken to metadata
    MultiAssayExperiment::metadata(expomicset)$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$steps,
      "run_pca"
    )

    return(expomicset)
  }else if (action=="get"){
    return(list(
      pca_df = tibble(dat),
      pca_feature = pca_feature,
      pca_sample = pca_sample,
      outliers = rownames(pca_sample$x)[outliers]
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
