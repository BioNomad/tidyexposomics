#' Plot Top Features by Factor from Integration Results
#'
#' Visualizes the top loading features for each factor from multi-omics integration results (e.g., MOFA or MCIA).
#'
#' @param expomicset A `MultiAssayExperiment` object containing integration results in the `metadata` slot (must include `integration_results`).
#' @param factors Character vector of factors to include (e.g., `"Factor1"`, `"Factor2"`). If `NULL`, all factors are plotted.
#' @param top_n Integer specifying the number of top features to show per factor. Default is `5`.
#' @param facet_cols Optional color palette for facet strip backgrounds (one per `exp_name`), used to distinguish factors.
#' @param exp_name_cols Optional color palette for experiment labels in the plot (`exp_name`), passed to `scale_color_manual()`.
#' @param alpha Numeric value between 0 and 1 controlling the transparency of facet strip background fill. Default is `0.5`.
#'
#' @details
#' This function supports the following integration methods:
#' - `"MOFA"`: Uses feature weights from MOFA2 (`get_weights()`).
#' - `"MCIA"`: Uses block loadings from MCIA (`@block_loadings`).
#'
#' For each factor, it:
#' - Selects the top `top_n` features by **absolute loading**.
#' - Creates a point-range plot showing the loading magnitude.
#' - Facets each factor with a customizable strip background.
#'
#' If palettes are not provided, defaults are chosen using `ggpubr::get_palette()`.
#'
#' @return A `ggplot2` object with one facet per factor, showing the top features and their loadings by experiment.
#'
#' @examples
#' \dontrun{
#' plot_top_factor_features(expomicset, factors = c("Factor1", "Factor2"), top_n = 10)
#' }
#'
#' @export
plot_top_factor_features <- function(
    expomicset,
    factors=NULL,
    top_n=5,
    facet_cols=NULL,
    exp_name_cols=NULL,
    alpha=0.5){

  require(ggplot2)

  # Check to see if multiomics integration results are available
  if(!"integration_results" %in% names(MultiAssayExperiment::metadata(expomicset)$multiomics_integration)){
    stop("Please run `multiomics_integration()` first.")
  }

  if(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method == "MOFA"){
    # Grab top MOFA features per factor
    df <- MOFA2::get_weights(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result) |>
      purrr::map2(
        names(MOFA2::get_weights(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result)),
        ~.x|>
          as.data.frame() |>
          tibble::rownames_to_column("feature") |>
          tidyr::pivot_longer(-feature,
                       names_to="factor",
                       values_to="loading") |>
          dplyr::mutate(abs_loading=abs(loading)) |>
          dplyr::mutate(exp_name=.y)) |>
      dplyr::bind_rows()

    # Filter by factors if provided
    if (!is.null(factors)) {
      df <- df |>
        dplyr::filter(factor %in% factors)
    }

    # Select top features per factor
    df <- df |>
      dplyr::group_by(factor) |>
      dplyr::arrange(dplyr::desc(abs_loading)) |>
      dplyr::slice_head(n=top_n) |>
      dplyr::mutate(feature=gsub("_.*","",feature))

    # If no facet_cols provided, use default
    if (!is.null(facet_cols)) {
      facet_cols <- facet_cols
    } else{
      facet_cols <- tidy_exp_pal[1:length(unique(df$factor))]
    }

    # If no exp_name_cols provided, use default
    if(!is.null(exp_name_cols)){
      exp_name_cols <- exp_name_cols
    } else{
      exp_name_cols <- rev(tidy_exp_pal)[1:length(unique(df$exp_name))]
    }

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
      theme(strip.text.y = element_text(face="bold.italic",angle = 0))+
      ggh4x::facet_grid2(
        factor~.,
        scales = "free_y",
        space = "free_y",
        strip = ggh4x::strip_themed(
          background_y = ggh4x::elem_list_rect(
            fill = scales::alpha(
              facet_cols,
              alpha
            )
          )
        )
      ) +
      scale_color_manual(values=exp_name_cols)+
      #facet_grid(factor~., scales="free_y")+
      #ggsci::scale_color_aaas()+
      labs(
        x="Absolute loading",
        y="",
        color="Experiment"
      )

  }else if(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method == "MCIA"){
    # Grab top MCIA features per factor
    df <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result@block_loadings |>
      purrr::map2(
        names(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result@block_loadings),
        ~.x|>
          as.data.frame() |>
          tibble::rownames_to_column("feature") |>
          tidyr::pivot_longer(-feature,
                       names_to="factor",
                       values_to="loading") |>
          dplyr::mutate(abs_loading=abs(loading)) |>
          dplyr::mutate(exp_name=.y)) |>
      dplyr::bind_rows()

    # Filter by factors if provided
    if (!is.null(factors)) {
      df <- df |>
        dplyr::filter(factor %in% factors)
    }

    # Select top features per factor
    df <- df |>
      dplyr::group_by(factor) |>
      dplyr::arrange(dplyr::desc(abs_loading)) |>
      dplyr::slice_head(n=top_n)

    # If no facet_cols provided, use default
    if (!is.null(facet_cols)) {
      facet_cols <- facet_cols
    } else{
      facet_cols <- tidy_exp_pal[1:length(unique(df$factor))]
    }

    # If no exp_name_cols provided, use default
    if(!is.null(exp_name_cols)){
      exp_name_cols <- exp_name_cols
    } else{
      exp_name_cols <- rev(tidy_exp_pal)[1:length(unique(df$exp_name))]
    }


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
      theme(strip.text.y = element_text(face="bold.italic",angle = 0))+
      ggh4x::facet_grid2(
        factor~.,
        scales = "free_y",
        space = "free_y",
        strip = ggh4x::strip_themed(
          background_y = ggh4x::elem_list_rect(
            fill = scales::alpha(
              facet_cols,
              alpha
            )
          )
        )
      ) +
      scale_color_manual(values=exp_name_cols)+
      # facet_grid(factor~., scales="free_y")+
      # ggsci::scale_color_npg()+
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
