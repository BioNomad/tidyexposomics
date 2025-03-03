plot_exposure_corr <- function(
    expomicset, 
    correlation_cutoff = 0.3) {
  require(igraph)
  require(ggraph)
  require(tidyverse)
  
  # Extract and filter relevant data
  correlation_data <- expomicset@metadata$exposure_correlation$filtered_table %>%
    filter(abs_correlation >= correlation_cutoff) 
  
  # Vertex Information
  vertex_df <- expomicset@metadata$var_info %>%
    filter(variable %in% c(correlation_data$var1, correlation_data$var2))
  
  # Create igraph object
  graph <- graph_from_data_frame(correlation_data,
                                 vertices = vertex_df,
                                 directed = FALSE)
  
  # Plot using ggraph
  ggraph(graph,
         layout = "linear", 
         circular = TRUE) + 
    geom_edge_arc(aes(
      color = correlation, 
      width = abs_correlation)) +  # Edge color & thickness
    geom_node_text(aes(
      label = name, 
      fontface = "bold.italic",
      angle = node_angle(x, y)),
      hjust = "outward",check_overlap = TRUE) +  # Node labels
    geom_node_point(shape = 21,
                    size = 4,
                    alpha=0.8,
                    aes(fill = category)) +  # Nodes
    theme_graph() +
    scale_fill_npg()+
    scale_edge_color_gradient2(
      low="blue", 
      mid="white", 
      high="red", 
      midpoint=0)+
    coord_fixed(xlim = c(-2, 2),
                ylim = c(-2, 2)) +
    #theme(text = element_text(face = "bold"))+
    guides(edge_width = "none")+
    labs(
      edge_color = "Correlation",
      fill = "Category"
    )
}





# a=graph_from_data_frame(expom@metadata$exposure_correlation$filtered_table,vertices = expom@metadata$var_info |> dplyr::select(variable,category))
# 
# ggraph(a, layout = "linear", circular = TRUE) + 
#   geom_edge_arc(aes(color = abs_correlation, width = abs_correlation)) +  # Edge color & thickness
#   geom_node_text(aes(label = name, angle = node_angle(x, y)), hjust = "outward") +  # Node labels
#   geom_node_point(shape = 21, size = 4, aes(fill = category)) +  # Nodes
#   theme_graph() +
#   scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0)+
#   coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
#   theme(axis.text = element_text(face = "bold"))+
#   guides()  # Hide node fill legend




# plot_exposure_corr <- function(expom, correlation_cutoff = 0.3) {
#   require(circlize)
#   require(dplyr)
#   require(viridis)
#   
#   # Extract and filter relevant data
#   correlation_data <- expom@metadata$exposure_correlation$filtered_table %>%
#     filter(abs_correlation >= correlation_cutoff) %>%
#     select(var1, var2, abs_correlation)
#   
#   # Define color mapping for absolute correlation values
#   col_fun <- colorRamp2(
#     c(min(correlation_data$abs_correlation), 
#       median(correlation_data$abs_correlation), 
#       max(correlation_data$abs_correlation)), 
#     c("blue","white","red")  # Smooth color gradient
#   )
#   
#   # Prepare data for chord diagram
#   mat <- correlation_data %>% 
#     spread(var2, abs_correlation, fill = 0) %>% 
#     column_to_rownames("var1") %>% 
#     as.matrix()
#   
#   # Define sector colors (optional customization)
#   sectors <- unique(c(correlation_data$var1, correlation_data$var2))
#   grid.col <- setNames(viridis(length(sectors)), sectors)
#   
#   # Clear previous Circos settings
#   circos.clear()
#   
#   # Generate chord diagram with proper track allocation
#   chordDiagram(
#     mat, 
#     grid.col = grid.col, 
#     col = col_fun(mat),  # Map colors to abs_correlation
#     annotationTrack = NULL, 
#     preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat)))))
#   )
#   
#   # Customize sector labels
#   circos.track(track.index = 1, panel.fun = function(x, y) {
#     circos.text(
#       CELL_META$xcenter, CELL_META$ylim[1] - 0.01,  # Move labels slightly outward
#       CELL_META$sector.index, 
#       facing = "clockwise", 
#       niceFacing = TRUE, 
#       adj = c(0, 0.5)
#     )
#   }, bg.border = NA)  # Remove background borders
# }