#' Plot Exposure-Outcome Associations
#'
#' Generates a forest plot visualizing significant exposure-outcome associations
#' from ExWAS results stored in `expomicset`.
#'
#' @param expomicset A `MultiAssayExperiment` object containing ExWAS results.
#' @param filter_col A character string specifying the column used for filtering significant associations. Default is `"p.value"`.
#' @param filter_thresh A numeric value specifying the threshold for filtering significant associations. Default is `0.05`.
#'
#' @details
#' The function extracts ExWAS results from `metadata(expomicset)$exwas$results_df`, 
#' filters based on `filter_col` and `filter_thresh`, and generates a forest plot of effect sizes.
#' Associations are color-coded as upregulated (red), downregulated (blue), or non-significant (grey).
#'
#' @return A `ggplot` object displaying the exposure-outcome associations.
#'
#' @examples
#' \dontrun{
#' plot_associate_exposure_outcome(expom)
#' }
#'
#' @export
plot_associate_exposure_outcome <- function(
    expomicset,
    filter_col = "p.value",
    filter_thresh = 0.05){
  require(ggplot2)
  
  # Check if the required metadata is present
  if(!"exwas" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `perform_exwas()` first.")
  }
  
  # Extract the results dataframe and filter based on the specified column and threshold
  exwas <- MultiAssayExperiment::metadata(expomicset)$exwas$results_df |> 
    dplyr::filter(!!sym(filter_col) < filter_thresh)
  
  covariates <- MultiAssayExperiment::metadata(expomicset)$exwas$covariates
  
  # Create forest plot of significant associations
  exwas  |> 
    dplyr::mutate(direction=case_when(
      estimate>0 & !!sym(filter_col) < filter_thresh ~ "up",
      estimate<0 & !!sym(filter_col) < filter_thresh ~ "down",
      .default = "ns"
    )) |> 
    ggplot(aes(x = estimate, 
               y = reorder(term,estimate),
               color = direction)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed") +
    geom_errorbarh(aes(
      xmin = estimate-std.error,
      xmax = estimate+std.error),
      color="grey55",
      height = 0.2) +
    geom_point(shape=18,
               size=5,
               alpha=0.5) +
    theme_bw() +
    scale_color_manual(values = c(
      "up" = "red3",
      "down" = "blue4",
      "ns" = "grey55"
    ))+
    theme(legend.position = "none",
          plot.subtitle = element_text(face="italic"))+
    labs(
      x = "Effect size",
      y = "",
      title = "Exposure-Outcome Associations",
      subtitle = paste("Covariates: ",paste(covariates, collapse = ", "))
      )
      
}

