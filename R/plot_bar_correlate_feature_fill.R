#' Plot Top Feature Associations Across Exposure Categories
#'
#' Generates a bar plot visualizing the most frequently associated features
#' across exposure categories based on correlation results.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use.
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#' @param top_n An integer specifying the number of top features to display. Default is `15`.
#'
#' @details
#' The function extracts correlation results from `metadata(expomicset)`,
#' groups features by exposure category, and visualizes the number of significant associations for each feature.
#' The top `top_n` features with the highest number of associations are displayed,
#' with color representing different exposure categories.
#'
#' @return A `ggplot` object displaying the top features and their associated exposure categories.
#'
#' @examples
#' \dontrun{
#' plot_bar_correlate_feature_fill(expom)
#' }
#'
#' @export
plot_bar_correlate_feature_fill <- function(
    expomicset,
    geneset = "degs",
    top_n = 15) {
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
    dplyr::group_by(feature,category) |>
    dplyr::reframe(n_category=length(category)) |>
    dplyr::group_by(feature) |>
    dplyr::mutate(total=sum(n_category)) |>
    dplyr::ungroup() |>
    dplyr::filter(feature %in% c(
      exp_feature_cor_df |>
        janitor::tabyl(feature) |>
        dplyr::arrange(desc(n)) |>
        dplyr::slice_head(n=15) |>
        dplyr::pull(feature)
    )) |>
    ggplot(aes(
      x=n_category,
      y=forcats::fct_reorder(feature, total),
      fill=category
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    scale_fill_tidy_exp()+
    scale_color_tidy_exp()+
    ggpubr::theme_pubr(legend="right",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"))+
    labs(title = "Frequently Associated Features",
         y = "",
         x = "No. of Associations",
         fill = "Exposure Category"
    )
}
