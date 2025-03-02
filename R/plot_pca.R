plot_pca <- function(
    expOmicSet,
    feature_col = "#00a9b2",
    sample_col = "#8a4f77",
    sample_outlier_col = "firebrick"
    
){
  require(tidyverse)
  require(ggfortify)
  require(ggpubr)
  require(ggsci)
  require(ggrepel)
  require(factoextra)
  require(patchwork)
  
  
  # grab data
  dat <- expOmicSet@metadata$pca$pca_df |> 
    as.data.frame()
  pca_feature <- expOmicSet@metadata$pca$pca_feature
  pca_sample <- expOmicSet@metadata$pca$pca_sample
  
  # capitalize exposure for plot
  dat <- dat |> 
    mutate(category = case_when(
    category == "exposure" ~ "Exposure",
    .default = category
  ))

  
  # PCA Feature Scatter Plot
  pca_plot_feature <- autoplot(pca_feature, 
                               data = dat, 
                               colour = "category") +
    theme_pubr(legend="right") +
    scale_color_npg() +
    labs(title = "PCA of Feature Space", color = "")
  
  # PCA Feature Scree Plot
  scree_feature <- fviz_eig(
    pca_feature,
    barfill = feature_col,
    barcolor = feature_col,
    main = "Scree Plot of Feature Space"
  )
  
  # Define outliers
  outlier_samples <- expOmicSet@metadata$pca$outliers
  
  # PCA sample scatter plot
  pca_plot_sample <- autoplot(pca_sample, colour = sample_col) +
    geom_text_repel(aes(
      x=PC1,
      y=PC2,
      label=ifelse(rownames %in% outlier_samples, rownames, NA)),
    color = sample_outlier_col) +
    theme_pubr() +
    labs(title = "PCA of Sample Space")
  
  
  # Sample Scree Plot
  scree_sample <- fviz_eig(
    pca_sample,
    barfill = sample_col,
    barcolor = sample_col,
    main = "Scree Plot of Sample Space"
  )
  
  # Combine plots
  combined_plot <- wrap_plots(
    pca_plot_feature, scree_feature,
    pca_plot_sample, scree_sample
  )
  
  return(combined_plot)

}
