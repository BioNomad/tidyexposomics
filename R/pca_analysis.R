pca_analysis <- function(expOmicSet) {
  require(tidyverse)
  require(ggfortify)
  require(ggpubr)
  require(ggsci)
  require(factoextra)
  require(patchwork)
  require(Hmisc)
  require(reshape2)
  
  # Identify common samples across all data
  message("Identifying common samples...")
  
  common_samples <- rownames(colData(expOmicSet))
  for (omics_name in names(experiments(expOmicSet))) {
    common_samples <- intersect(common_samples, colnames(experiments(expOmicSet)[[omics_name]]))
  }
  
  if (length(common_samples) == 0) {
    stop("No common samples found across exposure and omics data.")
  }
  
  # Subset colData to common samples
  message("Subsetting exposure data...")
  
  exposure_data <- colData(expOmicSet)[common_samples, ] |>
    as.data.frame() |>
    dplyr::select(where(is.numeric)) |>
    t() |>
    as.data.frame()
  exposure_data <- transform(exposure_data, category = "exposure")
  
  # Subset omics data to common samples
  message("Subsetting omics data...")
  
  omics_data <- lapply(names(experiments(expOmicSet)), function(omics_name) {
    assays(experiments(expOmicSet)[[omics_name]])[[1]][, common_samples, drop = FALSE] |>
      as.data.frame() |>
      transform(category = omics_name)
  }) |>
    bind_rows()
  
  # Combine datasets
  dat <- bind_rows(exposure_data, omics_data)
  dat <- transform(dat, id = rownames(dat))
  
  # Remove zero-variance columns
  feature_data <- dat |>
    dplyr::select(-c(category, id)) |>
    dplyr::select(where(\(x) var(x, na.rm = TRUE) > 0)) |> 
    mutate(across(where(is.numeric), ~ log2(.+abs(min(.))+1)))
  
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
  sample_data <- dplyr::select(sample_data, names(sample_vars[sample_vars > 0])) |> 
    mutate(across(where(is.numeric), ~ log2(.+abs(min(.))+1)))
  
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
  
  # Store results
  metadata(expOmicSet)$pca <- list(
    pca_df = tibble(dat),
    pca_feature = pca_feature,
    pca_sample = pca_sample,
    outliers = rownames(pca_sample$x)[outliers]
  )
  
  return(expOmicSet)
}