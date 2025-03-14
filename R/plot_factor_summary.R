#' Plot Summary of Factor Contributions from Multi-Omics Integration
#'
#' Generates a summary plot of factor contributions from multi-omics integration results.
#'
#' @param expomicset A `MultiAssayExperiment` object containing multi-omics integration results.
#'
#' @details
#' This function extracts integration results from `metadata(expomicset)$integration_results` and 
#' generates a summary plot of factor contributions. The visualization method depends on the 
#' integration approach used:
#'
#' - **MOFA**: Uses `MOFA2::plot_variance_explained()` to display variance explained by factors.
#' - **MCIA**: Uses `nipalsMCIA::block_weights_heatmap()` to show block loadings.
#'
#' The function automatically selects the appropriate visualization based on the integration method.
#'
#' @return A `ggplot` object displaying factor contributions for MOFA or a block weight heatmap for MCIA.
#'
#' @examples
#' \dontrun{
#' plot_factor_summary(expom)
#' }
#'
#' @export
plot_factor_summary <- function(
    expomicset
    ){
  if(!"integration_results" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `multiomics_integration()` first.")
  }
  
  if(MultiAssayExperiment::metadata(expomicset)$integration_results$method == "MOFA"){
    factor_contrib_plot <- MOFA2::plot_variance_explained(
      MultiAssayExperiment::metadata(expomicset)$integration_results$result,
      x="view",
      y="factor")+
      ggpubr::rotate_x_text(angle = 45)
    
  }else if(MultiAssayExperiment::metadata(expomicset)$integration_results$method == "MCIA"){
    factor_contrib_plot <- nipalsMCIA::block_weights_heatmap(
      MultiAssayExperiment::metadata(expomicset)$integration_results$result)
  }else{
    stop("Method not supported.")
  }
  return(factor_contrib_plot)
  
}
