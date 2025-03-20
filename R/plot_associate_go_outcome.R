#' Plot GO Group-Outcome Associations
#'
#' Generates a forest plot visualizing significant associations between GO groups and outcomes
#' from principal component regression results stored in `expomicset`.
#'
#' @param expomicset A `MultiAssayExperiment` object containing `pc_glm_results`.
#' @param filter_col A character string specifying the column used for filtering significant associations. Default is `"p.value"`.
#' @param filter_thresh A numeric value specifying the threshold for filtering significant associations. Default is `0.05`.
#' @param direction_filter A character string specifying which effect direction to display.
#' Options are `"all"`, `"up"` (positive associations), or `"down"` (negative associations). Default is `"all"`.
#' @param nrow An integer specifying the number of rows in the faceted plot layout. Default is `1`.
#'
#' @details
#' The function extracts GO group-outcome associations from `metadata(expomicset)$pc_glm_results`,
#' filters based on `filter_col` and `filter_thresh`, and visualizes the effect sizes using a forest plot.
#' Associations are color-coded as upregulated (red), downregulated (blue), or non-significant (grey).
#' Faceting is applied to separate results by experiment and cluster.
#'
#' @return A `ggplot` object displaying the GO group-outcome associations.
#'
#' @examples
#' \dontrun{
#' plot_associate_go_outcome(expom)
#' }
#'
#' @export
plot_associate_go_outcome <- function(
    expomicset,
    filter_col = "p.value",
    filter_thresh = 0.05,
    direction_filter = "all",
    nrow = 1){
  require(ggplot2)

  # Check if pc_glm_results exists in metadata
  if(!"pc_glm_results" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run the necessary analysis to generate `pc_glm_results` first.")
  }

  # Extract results
  results_df <- MultiAssayExperiment::metadata(expomicset)$pc_glm_results

  # Filter for significant associations
  filtered_df <- results_df |>
    dplyr::filter(!!sym(filter_col) < filter_thresh)

  # Create direction variable
  filtered_df <- filtered_df |>
    dplyr::mutate(direction=case_when(
      estimate > 0 & !!sym(filter_col) < filter_thresh ~ "up",
      estimate < 0 & !!sym(filter_col) < filter_thresh ~ "down",
      .default = "ns"
    ))

  # Apply direction filter
  if(direction_filter == "up"){
    filtered_df <- filtered_df |>
      dplyr::filter(direction == "up")
  } else if(direction_filter == "down"){
    filtered_df <- filtered_df |>
      dplyr::filter(direction == "down")
  }

  # Generate forest plot
  filtered_df |>
    tidyr::separate(col = "term",
             into=c("PC","exp_name","Cluster","go_group"),
             sep = "/") |>
    ggplot( aes(x = estimate,
                y = reorder(go_group, estimate),
                color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = estimate - std.error,
                       xmax = estimate + std.error),
                   color = "grey55", height = 0.2) +
    geom_point(shape = 18, size = 5, alpha = 0.5) +
    theme_bw() +
    scale_color_manual(values = c(
      "up" = "#8E0152",
      "down" = "#006666",
      "ns" = "grey55"
    )) +
    theme(legend.position = "none",
          plot.subtitle = element_text(face = "italic")) +
    ggh4x::facet_nested_wrap(exp_name+Cluster~.,nrow = nrow,scales = "free_x")+
    theme(strip.text = element_text(face="bold.italic"))+
    labs(
      x = "Effect size",
      y = "",
      title = "Go Group - Outcome Associations"
    )
}

