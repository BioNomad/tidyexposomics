#' Plot Exposure Associations by Category
#'
#' Generates a bar plot visualizing the number of significant exposure associations
#' across different exposure categories.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use.
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#'
#' @details
#' The function extracts correlation results from `metadata(expomicset)`,
#' groups exposures by category, and visualizes the number of significant associations.
#' The color represents different exposure categories.
#'
#' @return A `ggplot` object displaying the number of significant exposure associations per category.
#'
#' @examples
#' \dontrun{
#' plot_bar_correlate_exposure(expom)
#' }
#'
#' @export
plot_bar_correlate_exposure <- function(
    expomicset,
    geneset = "degs") {
  require(ggplot2)

  if(geneset=="degs"){
    # Check if the metadata contains the required data
    if(!"omics_exposure_deg_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation

  }else if(geneset=="factors"){
    # Check if the metadata contains the required data
    if(!"omics_exposure_factor_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }

  # Plot the bar plot
  exp_feature_cor_df |>
    janitor::tabyl(category) |>
    ggplot(aes(
      x=n,
      y=reorder(category,n),
      fill=category
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    geom_segment(aes(
      x = n,
      xend = n,
      y = as.numeric(reorder(category, n)) - 0.45,
      yend = as.numeric(reorder(category, n)) + 0.45,
      color = category,
    ), size = 1) +
    scale_fill_tidy_exp()+
    scale_color_tidy_exp()+
    ggpubr::theme_pubr(legend="none",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"))+
    labs(
      y = "",
      x = "No. of Associations",
      fill = "Exposure Category",
      title = "Exposure Associations")
}
