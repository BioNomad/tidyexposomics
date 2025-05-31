#' Plot Omics Feature Associations
#'
#' Generates a bar plot visualizing the number of significant associations between
#' omics features and exposures across different assays.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use.
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#'
#' @details
#' The function extracts correlation results from `metadata(expomicset)`,
#' groups features by assay, and visualizes the number of significant associations.
#' The color represents different assay types.
#'
#' @return A `ggplot` object displaying the number of significant feature associations per assay.
#'
#' @examples
#' \dontrun{
#' plot_bar_correlate_feature(expom)
#' }
#'
#' @export
plot_bar_correlate_feature <- function(
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

  exp_feature_cor_df |>
    janitor::tabyl(exp_name) |>
    ggplot(aes(
      x=n,
      y=reorder(exp_name,n),
      fill=exp_name
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    geom_segment(aes(
      x = n,
      xend = n,
      y = as.numeric(reorder(exp_name, n)) - 0.45,
      yend = as.numeric(reorder(exp_name, n)) + 0.45,
      color = exp_name,
    ), size = 1) +
    scale_fill_tidy_exp(rev = T)+
    scale_color_tidy_exp(rev = T)+
    ggpubr::theme_pubr(legend="none",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"))+
    labs(
      y = "",
      x = "No. of Associations",
      fill = "Exposure Category",
      title = "Omics Associations")
}
