plot_go_group_pc_outcome_assoc <- function(
    expomicset,
    filter_col = "p.value",
    filter_thresh = 0.05,
    direction_filter = "all"){  
  
  # Check if pc_glm_results exists in metadata
  if(!"pc_glm_results" %in% names(expomicset@metadata)){
    stop("Please run the necessary analysis to generate `pc_glm_results` first.")
  }
  
  # Extract results
  results_df <- expomicset@metadata$pc_glm_results
  
  # Filter for significant associations
  filtered_df <- results_df |> 
    filter(!!sym(filter_col) < filter_thresh)
  
  # Create direction variable
  filtered_df <- filtered_df |> 
    mutate(direction=case_when(
      estimate > 0 & !!sym(filter_col) < filter_thresh ~ "up",
      estimate < 0 & !!sym(filter_col) < filter_thresh ~ "down",
      .default = "ns"
    ))
  
  # Apply direction filter
  if(direction_filter == "up"){
    filtered_df <- filtered_df |> filter(direction == "up")
  } else if(direction_filter == "down"){
    filtered_df <- filtered_df |> filter(direction == "down")
  }
  
  # Generate forest plot
  ggplot(filtered_df, aes(x = estimate, 
                          y = reorder(term, estimate),
                          color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = estimate - std.error,
                       xmax = estimate + std.error),
                   color = "grey55", height = 0.2) +
    geom_point(shape = 18, size = 5, alpha = 0.5) +
    theme_bw() +
    scale_color_manual(values = c(
      "up" = "red3",
      "down" = "blue4",
      "ns" = "grey55"
    )) +
    theme(legend.position = "none",
          plot.subtitle = element_text(face = "italic")) +
    labs(
      x = "Effect size",
      y = "",
      title = "Exposure-Outcome Associations"
    )
}

