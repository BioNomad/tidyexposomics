#' Plot Summary of Exposure-Feature Correlations
#'
#' Generates a summary of exposure-feature associations using multiple bar plots,
#' combining feature and exposure-level insights.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use. 
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#' @param top_n An integer specifying the number of top exposures and features to display in the summary. Default is `15`.
#'
#' @details
#' This function compiles four bar plots:
#' - `plot_bar_correlate_exposure_fill()`: Top exposures associated with features.
#' - `plot_bar_correlate_feature_fill()`: Top features associated with exposures.
#' - `plot_bar_correlate_exposure()`: Summary of exposure category associations.
#' - `plot_bar_correlate_feature()`: Summary of omics assay associations.
#' 
#' The plots are arranged using `patchwork` to provide a high-level overview of exposure-feature relationships.
#'
#' @return A `patchwork` plot layout containing multiple bar plots summarizing exposure-feature associations.
#'
#' @examples
#' \dontrun{
#' plot_bar_correlate_summary(expom)
#' }
#'
#' @export
plot_bar_correlate_summary <- function(
    expomicset,
    geneset = "degs",
    top_n = 15) {
  require(ggplot2)
  require(patchwork)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  pop_exposures <- plot_bar_correlate_exposure_fill(expomicset,
                                      geneset = geneset,
                                      top_n)
  pop_features <- plot_bar_correlate_feature_fill(expomicset,
                                    geneset = geneset,
                                    top_n)
  exp_assoc <- plot_bar_correlate_exposure(expomicset,
                              geneset = geneset)
  omic_assoc <- plot_bar_correlate_feature(expomicset,
                                geneset = geneset)
  
  # Combine using patchwork
  ((pop_features/exp_assoc)+plot_layout(heights = c(3,1)) )| ((pop_exposures/omic_assoc)+patchwork::plot_layout(heights = c(3,1)))
}




