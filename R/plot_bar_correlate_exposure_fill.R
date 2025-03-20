#' Plot Top Exposure Associations Across Assays
#'
#' Generates a bar plot visualizing the most frequently associated exposures
#' across assays based on correlation results.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure-feature correlation results.
#' @param geneset A character string indicating the type of features to use.
#' Options are `"degs"` (differentially expressed genes) or `"factors"` (latent factors). Default is `"degs"`.
#' @param top_n An integer specifying the number of top exposures to display. Default is `15`.
#'
#' @details
#' The function extracts correlation results from `metadata(expomicset)`,
#' groups exposures by assay, and visualizes the number of significant associations for each exposure.
#' The top `top_n` exposures with the highest number of associations are displayed,
#' with color representing different assay types.
#'
#' @return A `ggplot` object displaying the top exposures and their associated assays.
#'
#' @examples
#' \dontrun{
#' plot_bar_correlate_exposure_fill(expom)
#' }
#'
#' @export
plot_bar_correlate_exposure_fill <- function(
    expomicset,
    geneset = "degs",
    top_n = 15) {
  require(ggplot2)

  if(geneset=="degs"){
    # Check if the metadata contains the required correlation data
    if(!"omics_exposure_deg_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_deg_correlation

  }else if(geneset=="factors"){
    # Check if the metadata contains the required correlation data
    if(!"omics_exposure_factor_correlation" %in% names(MultiAssayExperiment::metadata(expomicset))){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- MultiAssayExperiment::metadata(expomicset)$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }

  # Plot the bar plot
  exp_feature_cor_df |>
    dplyr::group_by(exposure,exp_name) |>
    dplyr::reframe(n_exp_name=length(exp_name)) |>
    dplyr::group_by(exposure) |>
    dplyr::mutate(total=sum(n_exp_name)) |>
    dplyr::ungroup() |>
    dplyr::filter(exposure %in% c(
      exp_feature_cor_df |>
        janitor::tabyl(exposure) |>
        dplyr::arrange(desc(n)) |>
        dplyr::slice_head(n=15) |>
        dplyr::pull(exposure)
    )) |>
    ggplot(aes(
      x=n_exp_name,
      y=forcats::fct_reorder(exposure, total),
      fill=exp_name
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    ggsci::scale_fill_aaas()+
    ggsci::scale_color_aaas()+
    ggpubr::theme_pubr(legend="right",base_size = 10)+
    theme(plot.title = element_text(face = "bold.italic"))+
    labs(title = "Frequently Associated Exposures",
         y = "",
         x = "No. of Associations",
         fill = "Assay Name"
    )
}
