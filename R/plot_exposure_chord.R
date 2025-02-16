
plot_exposure_chord <- function(expom, correlation_cutoff = 0.3) {
  require(circlize)
  require(dplyr)
  # Extract and filter relevant data
  correlation_data <- expom@metadata$exposure_correlation$filtered_table %>%
    filter(abs_correlation >= correlation_cutoff) %>%
    select(var1, var2, abs_correlation)
  
  # Define color mapping for absolute correlation values
  col_fun <- colorRamp2(
    c(min(correlation_data$abs_correlation), 
      median(correlation_data$abs_correlation), 
      max(correlation_data$abs_correlation)), 
    c("blue","white","red")  # Smooth color gradient
  )
  
  # Prepare data for chord diagram
  mat <- correlation_data %>% 
    spread(var2, abs_correlation, fill = 0) %>% 
    column_to_rownames("var1") %>% 
    as.matrix()
  
  # Define sector colors (optional customization)
  sectors <- unique(c(correlation_data$var1, correlation_data$var2))
  grid.col <- setNames(viridis(length(sectors)), sectors)
  
  # Clear previous Circos settings
  circos.clear()
  
  # Generate chord diagram with proper track allocation
  chordDiagram(
    mat, 
    grid.col = grid.col, 
    col = col_fun(mat),  # Map colors to abs_correlation
    annotationTrack = NULL, 
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat)))))
  )
  
  # Customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1] - 0.01,  # Move labels slightly outward
      CELL_META$sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE, 
      adj = c(0, 0.5)
    )
  }, bg.border = NA)  # Remove background borders
}

# Example Usage
plot_exposure_chord(expom, correlation_cutoff = 0.4)