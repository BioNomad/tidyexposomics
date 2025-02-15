plot_pca <- function(
    expOmicSet,
    feature_col = "#00a9b2",
    sample_col = "#8a4f77"
){
  require(tidyverse)
  require(ggpubr)
  
  return(expOmicSet@metadata$pca$combined_plot)

}