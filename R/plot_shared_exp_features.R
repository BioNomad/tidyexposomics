plot_shared_exp_features <- function(
    expomicset, 
    geneset = "degs",
    cutoff = 10) {
  require(igraph)
  require(ggraph)
  require(tidyverse)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expomicset@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expomicset@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expomicset@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expomicset@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  df <- expomicset@metadata$omics_exposure_deg_correlation
  
  exp_feature_cor_lst <- split(
    exp_feature_cor_df, 
    exp_feature_cor_df$exposure) |> 
    map(~ .x |>
          pull(feature) |>
          unique()) 
  
  exp_shared_feature_df <- .get_pairwise_overlaps(exp_feature_cor_lst) |> 
    filter(num_shared > cutoff) |> 
    arrange(num_shared)
  
  # Vertex Information
  vertex_df <- expomicset@metadata$var_info %>%
    filter(variable %in% c(
      exp_shared_feature_df$source,
      exp_shared_feature_df$target))
  
  # Create igraph object
  graph <- graph_from_data_frame(exp_shared_feature_df,
                                 vertices = vertex_df,
                                 directed = FALSE)
  
  # Plot using ggraph
  ggraph(graph,
         layout = "linear", 
         circular = TRUE) + 
    geom_edge_arc(aes(
      color = num_shared, 
      width = num_shared)) +  # Edge color & thickness
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
      midpoint = mean(exp_shared_feature_df$num_shared))+
    coord_fixed(xlim = c(-2, 2),
                ylim = c(-2, 2)) +
    #theme(text = element_text(face = "bold"))+
    guides(edge_width = "none")+
    labs(
      edge_width="Number of Shared Features",
      edge_color="Number of Shared Features",
      fill="Category"
    )
}

# 
# df <- expom@metadata$omics_exposure_deg_correlation
# 
# x <- split(df, df$exposure) |> 
#   map(~ .x |> pull(feature) |> unique()) 
# 
# y <- .get_pairwise_overlaps(x) |> 
#   filter(num_shared != 0)
