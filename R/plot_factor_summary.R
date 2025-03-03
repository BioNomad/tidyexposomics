plot_factor_summary <- function(
    expomicset
    ){
  require(MOFA2)
  require(nipalsMCIA)
  require(tidyverse)
  require(ggpubr)
  require(ggsci)
  if(!"integration_results" %in% names(expomicset@metadata)){
    stop("Please run `multiomics_integration()` first.")
  }
  
  if(expomicset@metadata$integration_results$method == "MOFA"){
    factor_contrib_plot <- plot_variance_explained(
      expomicset@metadata$integration_results$result,
      x="view",
      y="factor")+
      rotate_x_text(angle = 45)
    
  }else if(expomicset@metadata$integration_results$method == "MCIA"){
    factor_contrib_plot <- block_weights_heatmap(
      expomicset@metadata$integration_results$result)
  }else{
    stop("Method not supported.")
  }
  return(factor_contrib_plot)
  
}




# factor_contrib_plot <- block_weights_heatmap(expom@metadata$integration_results$result)
# 
# x <- expom@metadata$integration_results$result@block_loadings |> 
#   map2(names(expom@metadata$integration_results$result@block_loadings),~.x|> 
#         as.data.frame() |>
#         rownames_to_column("feature") |> 
#         pivot_longer(-feature, names_to="factor", values_to="loading") |> 
#         mutate(abs_loading=abs(loading)) |> 
#         mutate(exp_name=.y)) |> 
#   bind_rows() |> 
#   group_by(factor) |>
#   arrange(desc(abs_loading)) |> 
#   slice_head(n=5)
# 
# features_per_factor_plot <- x |> 
#   ggplot(aes(
#     x=abs_loading,
#     y=reorder(feature,abs_loading),
#     color=exp_name
#   ))+
#   geom_point(shape=18,
#              size=5,
#              alpha=0.5) +
#   geom_segment(aes(x=0, 
#                    xend=abs_loading, 
#                    y=feature, 
#                    yend=feature),
#                color="grey55")+
#   theme_bw()+
#   theme(strip.text = element_text(face="bold.italic"))+
#   facet_grid(factor~., scales="free_y")+
#   scale_color_npg()+
#   labs(
#     x="Absolute loading",
#     y="",
#     color="Experiment"
#   )
# 
# wrap_plots(features_per_factor_plot,
#            wrap_plots(grid.grabExpr(
#              draw(factor_contrib_plot)),
#            plot_spacer(),
#            plot_spacer(),plot_spacer(),
#            plot_spacer(),ncol=1))+
#   plot_layout(widths = c(1,3))
# 
# factor_contrib_plot <- plot_variance_explained(a@metadata$integration_results$result, x="view", y="factor")+
#   rotate_x_text(angle = 45)
# 
# y <- MOFA2::get_weights(a@metadata$integration_results$result) |> 
#   map2(names(MOFA2::get_weights(a@metadata$integration_results$result)),~.x|> 
#          as.data.frame() |>
#          rownames_to_column("feature") |> 
#          pivot_longer(-feature, names_to="factor", values_to="loading") |> 
#          mutate(abs_loading=abs(loading)) |> 
#          mutate(exp_name=.y)) |> 
#   bind_rows() |> 
#   group_by(factor) |>
#   arrange(desc(abs_loading)) |> 
#   slice_head(n=5) |> 
#   mutate(feature=gsub("_.*","",feature))
# 
# features_per_factor_plot <- y |> 
#   ggplot(aes(
#     x=abs_loading,
#     y=reorder(feature,abs_loading),
#     color=exp_name
#   ))+
#   geom_point(shape=18,
#              size=5,
#              alpha=0.5) +
#   geom_segment(aes(x=0, 
#                    xend=abs_loading, 
#                    y=feature, 
#                    yend=feature),
#                color="grey55")+
#   theme_bw()+
#   theme(strip.text = element_text(face="bold.italic"))+
#   facet_grid(factor~., scales="free_y")+
#   scale_color_npg()+
#   labs(
#     x="Absolute loading",
#     y="",
#     color="Experiment"
#   )
# 
# wrap_plots(features_per_factor_plot,
#            wrap_plots(factor_contrib_plot,
#              plot_spacer(),ncol=1))+
#   plot_layout(widths = c(1,3))
