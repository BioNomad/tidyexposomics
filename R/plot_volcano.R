#' Volcano Plot of Differential Abundance
#'
#' Generates a **volcano plot** to visualize differentially abundant features across assays.
#'
#' @param expomicset A `MultiAssayExperiment` object with differential abundance results.
#' @param pval_col A character string specifying the column for p-values. Default is `"adj.P.Val"`.
#' @param pval_thresh A numeric value for the significance threshold. Default is `0.05`.
#' @param logFC_col A character string specifying the column for log fold change. Default is `"logFC"`.
#' @param logFC_thresh A numeric value for the log fold change threshold. Default is `log2(1.5)`.
#' @param plot_n_sig A logical indicating whether to display the number of significant features in facet labels. Default is `TRUE`.
#' @param xlab Label for the x-axis. Default is `expression(Log[2]*"FC")`.
#' @param ylab Label for the y-axis. Default is `expression(-Log[10]*"P")`.
#' @param title Title of the plot. Default is `"Volcano Plot of Differential Abundance"`.
#' @param nrow Number of rows in the facet layout. Default is `2`.
#'
#' @details
#' This function:
#' - Extracts **differential abundance results** from `metadata(expomicset)`.
#' - Assigns significance based on `pval_thresh` and `logFC_thresh`.
#' - Colors points as **Upregulated (red)**, **Downregulated (blue)**, or **Not-Significant (gray)**.
#' - Adds **dashed cutoff lines** for significance thresholds.
#' - Facets by assay (`exp_name`), optionally displaying counts of significant features.
#'
#' @return A `ggplot2` object displaying a volcano plot of differential abundance.
#'
#' @examples
#' \dontrun{
#' plot_volcano(expom, pval_thresh = 0.01, logFC_thresh = log2(2))
#' }
#'
#' @export
plot_volcano <- function(
    expomicset,
    pval_col = "adj.P.Val",
    pval_thresh = 0.05,
    logFC_col = "logFC",
    logFC_thresh = log2(1.5),
    plot_n_sig = TRUE,
    xlab = expression(Log[2]*"FC"),
    ylab = expression(-Log[10]*"P"),
    title = "Volcano Plot of Differential Abundance",
    nrow=2
){
  require(ggplot2)
  
  # Check to see if Differential Abundance Results are available
  if(!"differential_abundance" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_differential_abundance()` first.")
  }
  
  if(plot_n_sig){
    # Grab the significant features per experiment
    exp_sum <- MultiAssayExperiment::metadata(expomicset)$differential_abundance |> 
      dplyr::group_by(exp_name) |> 
      dplyr::summarise(total = dplyr::n(),
                total_significant = sum(
                  !!sym(pval_col) < pval_thresh & 
                    abs(!!sym(logFC_col)) > logFC_thresh)) |> 
      dplyr::mutate(exp_name_plot=paste(
        exp_name,
        "\n",
        " (",
        total_significant,
        "/",
        total,
        ")",
        sep=""))
  } else {
    exp_sum <- MultiAssayExperiment::metadata(expomicset)$differential_abundance |> 
      dplyr::group_by(exp_name) |> 
      dplyr::summarise(total = dplyr::n(),
                total_significant = sum(
                  !!sym(pval_col) < pval_thresh & 
                    abs(!!sym(logFC_col)) > logFC_thresh)) |> 
      dplyr::mutate(exp_name_plot=exp_name)
  }
    
  
  MultiAssayExperiment::metadata(expomicset)$differential_abundance |> 
    dplyr::inner_join(exp_sum,
               by = "exp_name") |>
    dplyr::arrange(dplyr::desc(total)) |> 
    dplyr::mutate(exp_name_plot=factor(
      exp_name_plot,
      levels=unique(exp_name_plot))) |>
    dplyr::mutate(direction=dplyr::case_when(
      !!sym(logFC_col) > logFC_thresh &
        !!sym(pval_col) <  pval_thresh ~ "Upregulated",
      
      !!sym(logFC_col) < -logFC_thresh &
        !!sym(pval_col) <  pval_thresh ~ "Downregulated",
      
      .default = "Not-Significant"
    )) |> 
    ggplot(aes(
      x = !!sym(logFC_col),
      y = -log10(!!sym(pval_col)),
      color = direction)) +
    geom_point(alpha=0.5)+
    geom_vline(xintercept = c(-logFC_thresh, logFC_thresh),
               linetype = "dashed",
               color="grey55")+
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed",
               color="grey55")+
    ggpubr::theme_pubr(legend = "bottom")+
    scale_color_manual(values = c(
      "Upregulated" = "red3",
      "Downregulated" = "blue4",
      "Not-Significant" = "grey55"
    ))+
    facet_wrap(. ~ exp_name_plot,nrow=nrow)+
    theme(strip.text = element_text(face = "bold.italic"),
          plot.title = element_text(face = "bold.italic"))+
    labs(
      x = xlab,
      y = ylab,
      title = title, 
      color = "Direction"
    )
}

