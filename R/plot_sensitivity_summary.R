#' Plot Sensitivity Analysis Summary
#'
#' Generates a ridge plot and bar chart summarizing feature stability scores across assays.
#'
#' @param expomicset A `MultiAssayExperiment` object containing sensitivity analysis results
#' in `metadata(expomicset)$sensitivity_analysis`.
#' @param stability_score_thresh A numeric threshold for stability scores. Default is `NULL`,
#' which uses the threshold stored in `metadata(expomicset)$sensitivity_analysis$score_thresh`.
#' @param title A character string specifying the title of the ridge plot. Default is
#' `"Distribution of Stability Scores"`.
#'
#' @details
#' This function:
#' - Extracts feature stability scores from `metadata(expomicset)$sensitivity_analysis$feature_stability`.
#' - Displays a **ridge plot** of stability score distributions per assay.
#' - Displays a **bar chart** of the number of features per assay.
#' - Prints the number of features with stability scores above the threshold.
#'
#' @return A `patchwork` object combining a ridge plot and a bar chart.
#'
#' @examples
#' \dontrun{
#' plot_sensitivity_summary(expom)
#' }
#'
#' @export
plot_sensitivity_summary <- function(
    expomicset,
    stability_score_thresh = NULL,
    title="Distribution of Stability Scores"){

  if(!"sensitivity_analysis" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_sensitivity_analysis` first.")
  }

  require(ggplot2)
  require(patchwork)

  # Check if sensitivity analysis results are available
  if(!"sensitivity_analysis" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_sensitivity_analysis()` first.")
  }

  # Get the sensitivity scores
  sensitivity_sum <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$feature_stability |>
    dplyr::group_by(exp_name) |>
    dplyr::summarise(n=dplyr::n()) |>
    dplyr::arrange(desc(n))

  # Create a bar chart of the number of features per assay
  sensitivity_bar <- sensitivity_sum |>
    ggplot(aes(
      x=n,
      y=forcats::fct_reorder(exp_name, n),
      fill=exp_name
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    geom_segment(aes(
      x = n,
      xend = n,
      y = as.numeric(forcats::fct_reorder(exp_name, n)) - 0.45,
      yend = as.numeric(forcats::fct_reorder(exp_name, n)) + 0.45,
      color = exp_name,
    ), size = 1) +
    ggsci::scale_fill_aaas()+
    ggsci::scale_color_aaas()+
    ggpubr::theme_pubr(legend="none")+
    theme(plot.title = element_text(face = "bold.italic"),
          axis.text.y = element_blank(),
          text = element_text(size=10))+
    labs(title = "",
         y = "",
         x = "No. of Features"
         )

  # Identify the stability score threshold
  if(is.null(stability_score_thresh)){
    stability_score_thresh <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$score_thresh
  } else{
    stability_score_thresh <- stability_score_thresh
  }

  # Create a ridge plot of stability score distributions per assay
  sensitivity_ridgeplot <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$feature_stability |>
    dplyr::left_join(sensitivity_sum, by="exp_name") |>
    ggplot(aes(
      x=stability_score,
      y=forcats::fct_reorder(exp_name,n),
      fill=exp_name))+
    ggridges::geom_density_ridges()+
    ggsci::scale_fill_aaas()+
    theme_minimal()+
    geom_vline(xintercept=stability_score_thresh,
               linetype="dashed",
               color="grey55")+
    theme(plot.title = element_text(face = "bold.italic"),
          legend.position = "none")+
    labs(fill="Assay",
         x="Stability Score",
         y="",
         title = title)

  # Print the number of features with stability scores above the threshold
  sensitivity_sum <- MultiAssayExperiment::metadata(expomicset)$sensitivity_analysis$feature_stability |>
    dplyr::group_by(exp_name) |>
    dplyr::summarise(
      n_above=length(stability_score[stability_score>stability_score_thresh]),
      total=dplyr::n()) |>
    dplyr::arrange(desc(total))

  message("Number of Features with Stability Score > ",stability_score_thresh,":")
  for(i in 1:nrow(sensitivity_sum)){
    message(
      sensitivity_sum$exp_name[i],
      ": ",
      sensitivity_sum$n_above[i],
      " / ",
      sensitivity_sum$total[i])
  }

  return((sensitivity_ridgeplot|sensitivity_bar)+plot_layout(widths = c(3,1)))
}
