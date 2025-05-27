#' Plot Significant Exposure-Outcome Associations
#'
#' Generates a **forest plot** of significant associations between exposures and an outcome, based on results from `perform_exwas()`.
#'
#' @param expomicset A `MultiAssayExperiment` object that includes results from `perform_exwas()` in its metadata.
#' @param filter_col A character string specifying the column to use for significance filtering. Default is `"p.value"`.
#' @param filter_thresh A numeric threshold to apply to `filter_col`. Only terms below this threshold will be plotted. Default is `0.05`.
#' @param subtitle Optional subtitle for the plot. If `NULL`, the covariates used in the model will be shown.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts exposure-outcome association results from `metadata(expomicset)$exwas$results_df`.
#'   \item Filters associations using `filter_col < filter_thresh`.
#'   \item Classifies direction of effect as `"up"`, `"down"`, or `"ns"` (non-significant).
#'   \item Visualizes point estimates and standard errors using a horizontal forest plot.
#' }
#'
#' @return A `ggplot2` object showing a forest plot of significant exposure-outcome associations.
#'
#' @examples
#' \dontrun{
#' plot_associate_exposure_outcome(expomicset, filter_thresh = 0.01)
#' }
#'
#' @export
plot_associate_exposure_outcome <- function(
    expomicset,
    filter_col = "p.value",
    filter_thresh = 0.05,
    subtitle = NULL){
  require(ggplot2)

  # Check if the required metadata is present
  if(!"exwas" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `perform_exwas()` first.")
  }

  # Extract the results dataframe and filter based on the specified column and threshold
  exwas <- MultiAssayExperiment::metadata(expomicset)$exwas$results_df |>
    dplyr::filter(!!sym(filter_col) < filter_thresh)

  covariates <- MultiAssayExperiment::metadata(expomicset)$exwas$covariates

  # if subtitle is NULL, create a default subtitle
  if(is.null(subtitle)){
    subtitle <- paste("Covariates: ",paste(covariates, collapse = ", "))
  } else {
    subtitle <- subtitle
  }
  # Create forest plot of significant associations
  exwas  |>
    dplyr::mutate(direction=dplyr::case_when(
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
      "up" = "#8E0152",
      "down" = "#006666",
      "ns" = "grey55"
    ))+
    theme(legend.position = "none",
          plot.subtitle = element_text(face="italic"))+
    labs(
      x = "Effect size",
      y = "",
      title = "Exposure-Outcome Associations",
      subtitle = subtitle
      )

}

