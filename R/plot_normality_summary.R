#' Plot Normality Summary of Exposure Variables
#'
#' Generates a bar plot summarizing the number of exposure variables that pass or fail normality
#' tests (e.g., Shapiro-Wilk) before or after transformation.
#'
#' @param expomicset A `MultiAssayExperiment` object with quality control metadata.
#' @param transformed Logical; if `TRUE`, use results after transformation. Default is `FALSE`.
#'
#' @return A `ggplot` object summarizing the number of exposures classified as normal or not normal.
#'
#' @details
#' This function assumes that `run_normality_check()` has been executed and that the results are
#' stored in `metadata(expomicset)$quality_control$normality`. If `transformed = TRUE`, the function will
#' instead plot the transformation summary stored in `metadata(expomicset)$quality_control$transformation$norm_summary`,
#' which is populated by `transform_exposure()`.
#'
#' The plot includes both bar heights and overlaid line segments to reinforce the counts.
#' Colors are drawn from the Lancet or UChicago palettes via `ggsci` and `ggpubr`.
#'
#' @examples
#' \dontrun{
#' # Plot original exposure normality results
#' plot_normality_summary(expomicset)
#'
#' # Plot post-transformation normality summary
#' plot_normality_summary(expomicset, transformed = TRUE)
#' }
#'
#' @export
plot_normality_summary <- function(
    expomicset,
    transformed=FALSE
){
  require(ggplot2)

  # Check if "normality" is a name in metadata
  if(!("normality" %in% names(MultiAssayExperiment::metadata(expomicset)$quality_control))) {
    stop("Please run `run_normality_check() first.`")
  }

  # Check if "transformation" is a name in metadata
  if(transformed){
    if(!("transformation" %in% names(MultiAssayExperiment::metadata(expomicset)$quality_control))) {
      stop("Please run `transform_exposure() first.`")
    }
  }

  # Plot normality results
  if(transformed){
    norm_plot <- MultiAssayExperiment::metadata(expomicset)$quality_control$transformation$norm_summary |>
      ggplot(aes(
      x = var,
      y = value,
      fill = var
    )) +
      geom_bar(stat = "identity", alpha = 0.5) +
      geom_segment(
        aes(x = as.numeric(as.factor(var)) - 0.45,
            xend = as.numeric(as.factor(var)) + 0.45,
            y = value,
            yend = value,
            color = var),
        size = 1
      ) +
      ggpubr::theme_pubr(legend = "right") +
      scale_fill_manual(values = ggpubr::get_palette("uchicago",k = 2)[c(2,1)])+
      scale_color_manual(values = ggpubr::get_palette("uchicago",k = 2)[c(2,1)])+
      guides(color=F)+
      theme(plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic"))+
      labs(
        x = "",
        y = "No. of Exposures",
        fill = "",
        title = "Normality of Exposure Variables",
        subtitle = "Shapiro-Wilk Test"
      )
  }else{
    norm_plot <- MultiAssayExperiment::metadata(expomicset)$quality_control$normality$norm_summary |>
      ggplot(aes(
        x = var,
        y = value,
        fill = var
      )) +
      geom_bar(stat = "identity", alpha = 0.5) +
      geom_segment(
        aes(x = as.numeric(as.factor(var)) - 0.45,
            xend = as.numeric(as.factor(var)) + 0.45,
            y = value,
            yend = value,
            color = var),
        size = 1
      ) +
      ggpubr::theme_pubr(legend = "right") +
      ggsci::scale_fill_lancet() +
      ggsci::scale_color_lancet(guide = FALSE) +
      theme(plot.title = element_text(face = "bold.italic"),
            plot.subtitle = element_text(face = "italic"))+
      labs(
        x = "",
        y = "No. of Exposures",
        fill = "",
        title = "Normality of Exposure Variables",
        subtitle = "Shapiro-Wilk Test"
      )
  }

  # # Extract normality data
  # norm_df <- MultiAssayExperiment::metadata(expomicset)$normality$norm_df
  #
  # # Create normality plot
  # norm_plot <- table(norm_df$p.value > 0.05) |>
  #   as.data.frame() |>
  #   dplyr::mutate(Var1 = dplyr::case_when(
  #     Var1 == "FALSE" ~ "Not Normal",
  #     Var1 == "TRUE" ~ "Normal"
  #   )) |>
  #   ggplot(aes(
  #     x = Var1,
  #     y = Freq,
  #     fill = Var1
  #   )) +
  #   geom_bar(stat = "identity", alpha = 0.5) +
  #   geom_segment(
  #     aes(x = as.numeric(as.factor(Var1)) - 0.45,
  #         xend = as.numeric(as.factor(Var1)) + 0.45,
  #         y = Freq,
  #         yend = Freq,
  #         color = Var1),
  #     size = 1
  #   ) +
  #   ggpubr::theme_pubr(legend = "right") +
  #   ggsci::scale_fill_lancet() +
  #   ggsci::scale_color_lancet(guide = FALSE) +
  #   labs(
  #     x = "",
  #     y = "No. of Exposures",
  #     fill = "",
  #     title = "Normality of Exposure Variables",
  #     subtitle = "Shapiro-Wilk Test"
  #   )
  #
  # transform_plot <- MultiAssayExperiment::metadata(expomicset)$transformation$norm_summary |>
  #   dplyr::mutate(max_normal=ifelse(normal==max(normal),1,0)) |>
  #   dplyr::select(transformation, normal,not_normal,max_normal) |>
  #   tidyr::pivot_longer(cols = -c(transformation,max_normal),
  #                names_to = "normal",
  #                values_to = "Freq") |>
  #   dplyr::arrange(desc(max_normal)) |>
  #   dplyr::mutate(transformation=case_when(
  #     transformation=="log_trans" ~ "Log2",
  #     transformation=="sqrt_trans" ~ "Square Root",
  #     transformation=="x_1_3_trans" ~ "X^(1/3)",
  #     transformation=="none_trans" ~ "No Transformation"
  #   )) |>
  #   dplyr::mutate(normal=case_when(
  #     normal=="normal" ~ "Normal",
  #     normal=="not_normal" ~ "Not Normal"
  #   )) |>
  #   dplyr::mutate(transformation=factor(transformation)) |>
  #   ggplot(aes(
  #     x = reorder(transformation,-max_normal),
  #     y = Freq,
  #     fill = normal
  #   )) +
  #   geom_bar(stat = "identity", alpha = 0.5) +
  #   ggpubr::theme_pubr(legend = "right") +
  #   ggpubr::rotate_x_text(angle=45)+
  #   ggsci::scale_fill_lancet() +
  #   ggsci::scale_color_lancet(guide = FALSE) +
  #   labs(
  #     x = "",
  #     y = "No. of Exposures",
  #     fill = "",
  #     title = "Normality of Exposure Variables",
  #     subtitle = "Shapiro-Wilk Test"
  #   )

  return(norm_plot)

}


