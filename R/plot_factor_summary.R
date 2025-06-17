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
  if(!"integration_results" %in% names(MultiAssayExperiment::metadata(expomicset)$multiomics_integration)){
    stop("Please run `multiomics_integration()` first.")
  }

  require(ggplot2)

  if(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method == "MOFA"){
    factor_contrib_plot <- MOFA2::plot_variance_explained(
      MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$result,
      x="view",
      y="factor")+
      scale_fill_gradient2(low="#006666",
                           mid="white",
                           high="#8E0152",
                           midpoint = 0.5)+
      ggpubr::rotate_x_text(angle = 45)

  }else if(MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method == "MCIA"){
    # factor_contrib_plot <- nipalsMCIA::block_weights_heatmap(
    #   MultiAssayExperiment::metadata(expomicset)$integration_results$result)

    factor_contrib_plot <- expomicset |>
      MultiAssayExperiment::metadata() |>
      purrr::pluck("multiomics_integration") |>
      purrr::pluck("integration_results") |>
      purrr::pluck("result") |>
      purrr::pluck("block_score_weights") |>
      as.data.frame() |>
      tibble::rownames_to_column("omic") |>
      tidyr::pivot_longer(!omic,names_to = "factor",values_to = "weight") |>
      dplyr::mutate(factor=gsub("V","",factor),
             factor=factor(as.numeric(factor), levels = sort(unique(as.numeric(factor))))) |>
      ggplot(aes(x = factor,
                 y = omic,
                 fill = weight)) +
      geom_tile() +
      ggpubr::theme_pubr(legend = "right")+
      scale_fill_gradient2(low="#006666",
                           mid="white",
                           high="#8E0152",
                           midpoint = 0.5)+
      labs(
        x="Factor",
        y="",
        fill="Weight"
      )
  }else if (MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results$method == "MCCA"){

    expomicset |>
      MultiAssayExperiment::metadata() |>
      purrr::pluck("multiomics_integration") |>
      purrr::pluck("integration_results") |>
      purrr::pluck("result") |>
      map(~ {.x |> colMeans()}) |>
      bind_rows() |>
      mutate_all(~ abs(.)) |>
      rownames_to_column("factor") |>
      pivot_longer(-factor,names_to = "omic",values_to = "weight") |>
      dplyr::mutate(factor=factor(
        as.numeric(factor),
        levels = sort(unique(as.numeric(factor))))) |>
      # group_by(factor,.drop = T) |>
      # mutate(perc=weight/sum(weight)) |>
      ggplot(aes(x = factor,
                 y = omic,
                 fill = weight)) +
      geom_tile() +
      ggpubr::theme_pubr(legend = "right")+
      scale_fill_gradient2(low="#006666",
                           mid="white",
                           high="#8E0152",
                           midpoint = 0.5)+
      labs(
        x="Factor",
        y="",
        fill="Weight"
      )

  } else {
    stop("Method not supported.")
  }
  return(factor_contrib_plot)

}
