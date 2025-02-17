plot_volcano <- function(
    expOmicSet,
    pval_col = "adj.P.Val",
    pval_thresh = 0.05,
    logFC_col = "logFC",
    logFC_thresh = log2(1.5),
    plot_n_sig = TRUE,
    xlab = expression(Log[2]*"FC"),
    ylab = expression(-Log[10]*"P"),
    title = "Volcano Plot of Differential Abundance"
){
  require(tidyverse)
  require(ggpubr)
  
  if(!"differential_abundance" %in% names(expOmicSet@metadata)){
    stop("Please run `run_differential_abundance()` first.")
  }
  
  if(plot_n_sig){
    exp_sum <- expOmicSet@metadata$differential_abundance |> 
      group_by(assay_name) |> 
      summarise(total = n(),
                total_significant = sum(
                  !!sym(pval_col) < pval_thresh & 
                    abs(!!sym(logFC_col)) > logFC_thresh)) |> 
      mutate(assay_name_plot=paste(
        assay_name,
        "\n",
        " (",
        total_significant,
        "/",
        total,
        ")",
        sep=""))
  } else {
    exp_sum <- expOmicSet@metadata$differential_abundance |> 
      group_by(assay_name) |> 
      summarise(total = n(),
                total_significant = sum(
                  !!sym(pval_col) < pval_thresh & 
                    abs(!!sym(logFC_col)) > logFC_thresh)) |> 
      mutate(assay_name_plot=assay_name)
  }
    
  
    expOmicSet@metadata$differential_abundance |> 
    inner_join(exp_sum,
               by = "assay_name") |>
    arrange(desc(total)) |> 
    mutate(assay_name_plot=factor(
      assay_name_plot,
      levels=unique(assay_name_plot))) |>
    mutate(direction=case_when(
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
    theme_pubr(legend = "bottom")+
    scale_color_manual(values = c(
      "Upregulated" = "red3",
      "Downregulated" = "blue4",
      "Not-Significant" = "grey55"
    ))+
    facet_grid(. ~ assay_name_plot)+
    theme(strip.text = element_text(face = "bold.italic"),
          plot.title = element_text(face = "bold.italic"))+
    labs(
      x = xlab,
      y = ylab,
      title = title, 
      color = "Direction"
    )
}

