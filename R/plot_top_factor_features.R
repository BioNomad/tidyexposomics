plot_top_factor_features <- function(
    expOmicSet,
    top_n=5){
  
  if(!"integration_results" %in% names(expOmicSet@metadata)){
    stop("Please run `multiomics_integration()` first.")
  }
  
  if(expOmicSet@metadata$integration_results$method == "MOFA"){
    df <- MOFA2::get_weights(expOmicSet@metadata$integration_results$result) |> 
      map2(
        names(MOFA2::get_weights(expOmicSet@metadata$integration_results$result)),
        ~.x|> 
          as.data.frame() |>
          rownames_to_column("feature") |> 
          pivot_longer(-feature,
                       names_to="factor",
                       values_to="loading") |> 
          mutate(abs_loading=abs(loading)) |> 
          mutate(exp_name=.y)) |> 
      bind_rows() |> 
      group_by(factor) |>
      arrange(desc(abs_loading)) |> 
      slice_head(n=5) |> 
      mutate(feature=gsub("_.*","",feature))
    
    features_per_factor_plot <- df |> 
      ggplot(aes(
        x=abs_loading,
        y=reorder(feature,abs_loading),
        color=exp_name
      ))+
      geom_point(shape=18,
                 size=5,
                 alpha=0.5) +
      geom_segment(aes(x=0, 
                       xend=abs_loading, 
                       y=feature, 
                       yend=feature),
                   color="grey55")+
      theme_bw()+
      theme(strip.text = element_text(face="bold.italic"))+
      facet_grid(factor~., scales="free_y")+
      scale_color_npg()+
      labs(
        x="Absolute loading",
        y="",
        color="Experiment"
      )
    
  }else if(expOmicSet@metadata$integration_results$method == "MCIA"){
    df <- expOmicSet@metadata$integration_results$result@block_loadings |> 
      map2(
        names(expOmicSet@metadata$integration_results$result@block_loadings),
        ~.x|> 
          as.data.frame() |>
          rownames_to_column("feature") |> 
          pivot_longer(-feature,
                       names_to="factor", 
                       values_to="loading") |> 
          mutate(abs_loading=abs(loading)) |> 
          mutate(exp_name=.y)) |> 
      bind_rows() |> 
      group_by(factor) |>
      arrange(desc(abs_loading)) |> 
      slice_head(n=top_n)
    
    features_per_factor_plot <- df |> 
      ggplot(aes(
        x=abs_loading,
        y=reorder(feature,abs_loading),
        color=exp_name
      ))+
      geom_point(shape=18,
                 size=5,
                 alpha=0.5) +
      geom_segment(aes(x=0, 
                       xend=abs_loading, 
                       y=feature, 
                       yend=feature),
                   color="grey55")+
      theme_bw()+
      theme(strip.text = element_text(face="bold.italic"))+
      facet_grid(factor~., scales="free_y")+
      scale_color_npg()+
      labs(
        x="Absolute loading",
        y="",
        color="Experiment"
      )
    
  }else{
    stop("Method not supported.")
  }
  return(features_per_factor_plot)
  
}