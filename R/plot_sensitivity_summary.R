plot_sensitivity_summary <- function(
    expomicset,
    title="Distribution of Stability Scores"){
  
  if(!"sensitivity_analysis" %in% names(expomicset@metadata)){
    stop("Please run `run_sensitivity_analysis` first.")
  }
  
  require(ggpubr)
  require(ggsci)
  require(ggridges)
  require(patchwork)
  
  if(!"sensitivity_analysis" %in% names(expomicset@metadata)){
    stop("Please run `run_sensitivity_analysis()` first.")
  }
  
  sensitivity_sum <- expomicset@metadata$sensitivity_analysis$feature_stability |>
    group_by(exp_name) |> 
    dplyr::summarise(n=n()) |> 
    arrange(desc(n)) 
  
  sensitivity_bar <- sensitivity_sum |> 
    ggplot(aes(
      x=n,
      y=fct_reorder(exp_name, n),
      fill=exp_name
    ))+
    geom_bar(stat = "identity",alpha=0.7) +
    geom_segment(aes(
      x = n,                    
      xend = n,                    
      y = as.numeric(fct_reorder(exp_name, n)) - 0.45,
      yend = as.numeric(fct_reorder(exp_name, n)) + 0.45,
      color = exp_name,
    ), size = 1) +
    scale_fill_npg()+
    scale_color_npg()+
    theme_pubr(legend="none")+
    theme(plot.title = element_text(face = "bold.italic"),
          axis.text.y = element_blank(),
          text = element_text(size=10))+
    labs(title = "",
         y = "",
         x = "No. of Features"
         #x = paste("No. of Differentially", "Abundant Features",sep="\n")
         ) 
  
  score_thresh <- expomicset@metadata$sensitivity_analysis$score_thresh
  
  sensitivity_ridgeplot <- expomicset@metadata$sensitivity_analysis$feature_stability |>
    left_join(sensitivity_sum, by="exp_name") |>
    ggplot(aes(
      x=stability_score,
      y=fct_reorder(exp_name,n),
      fill=exp_name))+
    geom_density_ridges()+
    scale_fill_npg()+
    theme_minimal()+
    geom_vline(xintercept=score_thresh,
               linetype="dashed",
               color="grey55")+
    theme(plot.title = element_text(face = "bold.italic"),
          legend.position = "none")+
    labs(fill="Assay",
         x="Stability Score",
         y="",
         title = title)
  
  
  return((sensitivity_ridgeplot|sensitivity_bar)+plot_layout(widths = c(3,1)))
}
