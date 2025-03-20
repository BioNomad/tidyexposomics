#' Plot Top Features for Each Factor
#'
#' Visualizes the most influential features for each latent factor based on multi-omics integration.
#'
#' @param expomicset A `MultiAssayExperiment` object containing factor loadings from MOFA or MCIA.
#' @param top_n An integer specifying the number of top features to display per factor. Default is `5`.
#'
#' @details
#' This function:
#' - Extracts factor loadings from MOFA (`get_weights()`) or MCIA (`@block_loadings`).
#' - Selects the top `top_n` features per factor based on absolute loading values.
#' - Generates a **facet grid scatter plot** of the top features colored by assay type.
#'
#' **Supported Integration Methods:**
#' - **MOFA**: Uses `get_weights()` to extract factor loadings.
#' - **MCIA**: Uses `@block_loadings` to extract block loadings.
#'
#' @return A `ggplot2` object displaying a scatter plot of the top factor-associated features.
#'
#' @examples
#' \dontrun{
#' plot_top_factor_features(expom, top_n = 10)
#' }
#'
#' @export
plot_top_factor_features <- function(
    expomicset,
    top_n=5){

  require(ggplot2)

  # Check to see if multiomics integration results are available
  if(!"integration_results" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `multiomics_integration()` first.")
  }

  if(MultiAssayExperiment::metadata(expomicset)$integration_results$method == "MOFA"){
    # Grab top MOFA features per factor
    df <- MOFA2::get_weights(MultiAssayExperiment::metadata(expomicset)$integration_results$result) |>
      purrr::map2(
        names(MOFA2::get_weights(MultiAssayExperiment::metadata(expomicset)$integration_results$result)),
        ~.x|>
          as.data.frame() |>
          tibble::rownames_to_column("feature") |>
          tidyr::pivot_longer(-feature,
                       names_to="factor",
                       values_to="loading") |>
          dplyr::mutate(abs_loading=abs(loading)) |>
          dplyr::mutate(exp_name=.y)) |>
      dplyr::bind_rows() |>
      dplyr::group_by(factor) |>
      dplyr::arrange(dplyr::desc(abs_loading)) |>
      dplyr::slice_head(n=5) |>
      dplyr::mutate(feature=gsub("_.*","",feature))

    # Create a plot of top features per factor
    features_per_factor_plot <- df |>
      ggplot(aes(
        x=abs_loading,
        y=reorder(feature,abs_loading),
        color=exp_name
      ))+
      geom_point(shape=18,
                 size=5,
                 alpha=0.5) +
      geom_segment(aes(x=0,
                       xend=abs_loading,
                       y=feature,
                       yend=feature),
                   color="grey55")+
      theme_bw()+
      theme(strip.text = element_text(face="bold.italic"))+
      facet_grid(factor~., scales="free_y")+
      ggsci::scale_color_aaas()+
      labs(
        x="Absolute loading",
        y="",
        color="Experiment"
      )

  }else if(MultiAssayExperiment::metadata(expomicset)$integration_results$method == "MCIA"){
    # Grab top MCIA features per factor
    df <- MultiAssayExperiment::metadata(expomicset)$integration_results$result@block_loadings |>
      purrr::map2(
        names(MultiAssayExperiment::metadata(expomicset)$integration_results$result@block_loadings),
        ~.x|>
          as.data.frame() |>
          tibble::rownames_to_column("feature") |>
          tidyr::pivot_longer(-feature,
                       names_to="factor",
                       values_to="loading") |>
          dplyr::mutate(abs_loading=abs(loading)) |>
          dplyr::mutate(exp_name=.y)) |>
      dplyr::bind_rows() |>
      dplyr::group_by(factor) |>
      dplyr::arrange(dplyr::desc(abs_loading)) |>
      dplyr::slice_head(n=top_n)

    # Create a plot of top features per factor
    features_per_factor_plot <- df |>
      ggplot(aes(
        x=abs_loading,
        y=reorder(feature,abs_loading),
        color=exp_name
      ))+
      geom_point(shape=18,
                 size=5,
                 alpha=0.5) +
      geom_segment(aes(x=0,
                       xend=abs_loading,
                       y=feature,
                       yend=feature),
                   color="grey55")+
      theme_bw()+
      theme(strip.text = element_text(face="bold.italic"))+
      facet_grid(factor~., scales="free_y")+
      ggsci::scale_color_npg()+
      labs(
        x="Absolute loading",
        y="",
        color="Experiment"
      )

  }else{
    stop("Method not supported.")
  }
  return(features_per_factor_plot)

}
