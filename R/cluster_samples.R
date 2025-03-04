cluster_samples <- function(expomicset, 
                            exposure_cols = NULL,  
                            dist_method = NULL, 
                            cluster_method = "ward.D", 
                            clustering_approach = "diana",
                            action = "add") {
  require(MultiAssayExperiment)
  require(tidyverse)
  require(ComplexHeatmap)
  require(cluster)
  require(vegan)  
  require(circlize)
  require(FactoMineR) 
  require(factoextra)
  require(dynamicTreeCut)
  require(densityClust)
  
  message("Starting clustering analysis...")
  
  # Validate exposure_cols input
  if (is.null(exposure_cols) || length(exposure_cols) == 0) {
    stop("Please specify `exposure_cols` as a character vector of `colData` column names.")
  }
  
  # Ensure selected exposure columns exist in colData
  available_variables <- intersect(
    exposure_cols,
    colnames(MultiAssayExperiment::colData(expomicset)))
  
  if (length(available_variables) == 0) {
    stop("None of the specified `exposure_cols` are available in `colData`. Clustering cannot proceed.")
  }
  
  # Extract numeric exposure data for clustering
  exposure_data <- MultiAssayExperiment::colData(expomicset) |> 
    as.data.frame() |> 
    dplyr::select(dplyr::all_of(available_variables))
  
  # Remove missing/infinite values
  exposure_data <- exposure_data |> 
    na.omit() |> 
    dplyr::mutate(
      dplyr::across(dplyr::everything(),
                    ~ ifelse(is.infinite(.),
                             NA,
                             .))) |> 
    na.omit()
  
  if (nrow(exposure_data) == 0) {
    stop("All selected exposure data contains missing or infinite values.")
  }
  
  # Determine appropriate distance metric
  if (is.null(dist_method)) {
    if (all(sapply(exposure_data, is.numeric))) {
      dist_method <- "euclidean"  # Continuous data
    } else {
      dist_method <- "gower"  # Mixed data types
    }
  }
  
  # Compute distance matrix
  sample_dist <- if (dist_method == "gower") {
    cluster::daisy(exposure_data, metric = "gower")
  } else {
    dist(exposure_data, method = dist_method)
  }
  
  # Function to determine optimal k based on the selected clustering approach
  determine_k <- function(dist_matrix, cluster_method) {
    if (clustering_approach == "diana") {
      
      # Determine optimal k using the height difference method
      sample_cluster <- cluster::diana(as.dist(dist_matrix))
      height_diffs <- diff(sample_cluster$height)
      cutoff_index <- which.max(height_diffs)
      return(length(sample_cluster$height) - cutoff_index)
      
    } else if (clustering_approach == "gap") {
      
      # Determine optimal k using the gap statistic
      gap_stat <- cluster::clusGap(as.matrix(dist_matrix), FUN = hcut, K.max = 20, B = 50)
      return(cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"]))
      
    } else if (clustering_approach == "elbow") {
      
      # Determine optimal k using the elbow method
      wss_plot <- factoextra::fviz_nbclust(as.matrix(dist_matrix), FUN = hcut, method = "wss")
      
      # Identify the first significant drop and ensure it's a number we can use
      k_optimal <- which.min(diff(diff(wss_plot$data$y))) + 1  
      if (is.na(k_optimal) || k_optimal < 2) k_optimal <- 3 
      return(k_optimal)
      
    } else if (clustering_approach == "dynamic") {
      
      # Determine optimal k using dynamic tree cut
      sample_cluster <- hclust(as.dist(dist_matrix), 
                               method = cluster_method)
      cut_clusters <- dynamicTreeCut::cutreeDynamic(
        dendro = as.dendrogram(sample_cluster),
        distM = as.matrix(dist_matrix), deepSplit = 2)
      return(length(unique(cut_clusters)))
      
    } else if (clustering_approach == "density") {
      
      # Determine optimal k using density-based clustering
      dclust <- densityClust::densityClust(dist_matrix, gaussian = TRUE)
      dclust <- findClusters(dclust, rho = 0.3, delta = 0.5)
      return(max(dclust$clusters))
      
    } else {
      stop("Invalid clustering approach selected.")
    }
  }
  
  # Determine optimal k clusters
  k_samples <- determine_k(sample_dist, cluster_method)
  
  # Perform hierarchical clustering
  sample_cluster <- hclust(as.dist(sample_dist), method = cluster_method)
  
  # Cut dendrogram using optimal k
  sample_groups <- cutree(sample_cluster, k = k_samples)
  
  message("Optimal number of clusters for samples: ", k_samples)
  
  # Generate heatmap
   heatmap <- ComplexHeatmap::Heatmap(
    matrix = t(scale(exposure_data)),  # Scale and transpose for clustering
    name = "Scaled Exposure",
    cluster_columns = sample_cluster,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      df=data.frame(
      Cluster = as.character(sample_groups),
      id = names(sample_groups)) |> 
        column_to_rownames("id")),
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8, angle = 45),
    heatmap_legend_param = list(title = "Exposure\nValue"),
    col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
  
  if (action == "add") {
    # Save clustering results in metadata
    MultiAssayExperiment::metadata(expomicset)$sample_clustering <- list(
      sample_cluster = sample_cluster,
      sample_groups = sample_groups,
      heatmap = heatmap
    )
    return(expomicset)
  } else if (action == "get") {
    return(list(
      sample_cluster = sample_cluster,
      sample_groups = sample_groups,
      heatmap = heatmap
    ))
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
}
