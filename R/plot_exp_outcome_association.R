plot_exp_outcome_association <- function(
    expomicset,
    filter_col = "p.value",
    filter_thresh = 0.05){
  
  if(!"exwas" %in% names(expomicset@metadata)){
    stop("Please run `perform_exwas()` first.")
  }
  
  exwas <- expomicset@metadata$exwas$results_df |> 
    filter(!!sym(filter_col) < filter_thresh)
  
  covariates <- expomicset@metadata$exwas$covariates
  
  # create forest plot of significant associations
  exwas  |> 
    mutate(direction=case_when(
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

