plot_sensitivity_summary <- function(
    expOmicSet,
    title="Distribution of Stability Scores"){
  
  if(!"sensitivity_analysis" %in% names(expOmicSet@metadata)){
    stop("Please run `run_sensitivity_analysis` first.")
  }
  
  score_thresh <- expOmicSet@metadata$sensitivity_analysis$score_thresh
  
  p <- expOmicSet@metadata$sensitivity_analysis$feature_stability |>
    ggplot(aes(
      x=stability_score,
      y=exp_name,
      fill=exp_name))+
    geom_density_ridges()+
    scale_fill_npg()+
    theme_minimal()+
    geom_vline(xintercept=score_thresh,
               linetype="dashed",
               color="grey55")+
    theme(plot.title = element_text(face = "bold.italic"))+
    labs(fill="Assay",
         x="Stability Score",
         y="",
         title = title)
  
  return(p)
}
