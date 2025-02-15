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
  pca_plot_feature <- autoplot(pca_feature, data = dat, colour = "category") +
    theme_pubr(legend="right") +
    scale_color_npg() +
    labs(title = "PCA of Feature Space", color = "")
  scree_feature <- fviz_eig(
    pca_feature,
    barfill = "#00a9b2",
    barcolor = "#00a9b2",
    main = "Scree Plot of Feature Space"
  )
  
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
  pca_plot_sample <- autoplot(pca_sample, colour = "#8a4f77") +
    theme_pubr() +
    labs(title = "PCA of Sample Space")
  scree_sample <- fviz_eig(
    pca_sample,
    barfill = "#8a4f77",
    barcolor = "#8a4f77",
    main = "Scree Plot of Sample Space"
  )
  
  # Combine plots
  combined_plot <- wrap_plots(
    pca_plot_feature, scree_feature,
    pca_plot_sample, scree_sample
  )
  
  # Store results
  metadata(expOmicSet)$pca <- list(
    pca_feature = pca_feature,
    pca_plot_feature = pca_plot_feature,
    scree_feature = scree_feature,
    pca_sample = pca_sample,
    pca_plot_sample = pca_plot_sample,
    scree_sample = scree_sample,
    combined_plot = combined_plot
  )
  
  return(expOmicSet)
}