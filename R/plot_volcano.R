#' Volcano Plot of Differential Abundance
#'
#' Generates a **volcano plot** to visualize differential abundance results across one or more omics layers.
#'
#' @param expomicset A `MultiAssayExperiment` object containing differential abundance results in `metadata(expomicset)$differential_abundance`.
#' @param pval_col A character string specifying the column containing p-values. Default is `"adj.P.Val"`.
#' @param pval_thresh A numeric threshold for significance. Features with p-values below this are considered significant. Default is `0.05`.
#' @param logFC_col A character string specifying the column for log fold changes. Default is `"logFC"`.
#' @param logFC_thresh A numeric threshold for absolute log fold change significance. Default is `log2(1.5)`.
#' @param plot_n_sig Logical; if `TRUE`, appends the number of significant features to facet titles. Default is `TRUE`.
#' @param top_n_label Optional integer. If provided, the top `n` most significant features per assay will be labeled on the plot.
#' @param features_to_label Optional character vector. Specific features to label regardless of significance.
#' @param feature_col A character string naming the feature ID column to use for labeling. Default is `"feature"`.
#' @param xlab Label for the x-axis. Default is `expression(Log[2]*"FC")`.
#' @param ylab Label for the y-axis. Default is `expression(-Log[10]*"P")`.
#' @param title Plot title. Default is `"Volcano Plot of Differential Abundance"`.
#' @param nrow Number of rows in the `facet_wrap()` layout. Default is `2`.
#'
#' @details
#' The function:
#' - Extracts differential abundance results from `metadata(expomicset)$differential_abundance`.
#' - Assigns each feature a direction of change: **Upregulated**, **Downregulated**, or **Not-Significant**.
#' - Uses `logFC_thresh` and `pval_thresh` to define thresholds.
#' - Adds dashed lines to indicate cutoffs for fold change and significance.
#' - Uses `facet_wrap()` to display each assay (`exp_name`) separately.
#' - Optionally labels the most significant features or user-defined ones.
#'
#' @return A `ggplot2` object representing the volcano plot.
#'
#' @examples
#' \dontrun{
#' plot_volcano(expom, pval_thresh = 0.01, logFC_thresh = log2(2))
#' plot_volcano(expom, features_to_label = c("TP53", "MYC"))
#' plot_volcano(expom, top_n_label = 5)
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
    top_n_label = NULL,
    features_to_label = NULL,
    feature_col = "feature",
    xlab = expression(Log[2]*"FC"),
    ylab = expression(-Log[10]*"P"),
    title = "Volcano Plot of Differential Abundance",
    nrow=2
){
  require(ggplot2)

  # Check to see if Differential Abundance Results are available
  if(!"differential_abundance" %in% names(MultiAssayExperiment::metadata(expomicset)$differential_analysis)){
    stop("Please run `run_differential_abundance()` first.")
  }

  if(plot_n_sig){
    # Grab the significant features per experiment
    exp_sum <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance |>
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
    exp_sum <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance |>
      dplyr::group_by(exp_name) |>
      dplyr::summarise(total = dplyr::n(),
                total_significant = sum(
                  !!sym(pval_col) < pval_thresh &
                    abs(!!sym(logFC_col)) > logFC_thresh)) |>
      dplyr::mutate(exp_name_plot=exp_name)
  }


  plot_df <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance |>
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
    ))

  volcano <- plot_df |>
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
      "Upregulated" = "#8E0152",
      "Downregulated" = "#006666",
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

  if(!is.null(top_n_label)){
    # Label the top n points
    volcano <- volcano +
      ggrepel::geom_label_repel(
        data = plot_df |>
          group_by(exp_name) |>
          dplyr::arrange(!!sym(pval_col)) |>
          dplyr::slice_head(n = top_n_label),
        aes(label = paste0("italic('", !!sym(feature_col), "')")),
        parse = TRUE,
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      )
      # ggrepel::geom_label_repel(
      #   data = plot_df |>
      #     group_by(exp_name) |>
      #     dplyr::arrange(!!sym(pval_col)) |>
      #     dplyr::slice_head(n = top_n_label),
      #   aes(label = !!sym(feature_col)),
      #   size = 3,
      #   max.overlaps = Inf,
      #   show.legend = FALSE
      # )
  }

  if(!is.null(features_to_label)) {
    # Label specific features
    volcano <- volcano +
      ggrepel::geom_label_repel(
        data = plot_df |>
          dplyr::filter(!!sym(feature_col) %in% features_to_label),
        aes(label = paste0("italic('", !!sym(feature_col), "')")),
        parse = TRUE,
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      )
      # ggrepel::geom_label_repel(
      #   data = plot_df |>
      #     dplyr::filter(!!sym(feature_col) %in% features_to_label),
      #   aes(label = !!sym(feature_col)),
      #   size = 3,
      #   max.overlaps = Inf,
      #   show.legend = FALSE
      # )
  }

  volcano
}

