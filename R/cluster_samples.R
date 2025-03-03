cluster_samples <- function(expomicset, 
                            exposure_cols = NULL,  # Character vector of `colData` columns to cluster
                            dist_method = NULL, 
                            cluster_method = "ward.D", 
                            clustering_approach = "diana",
                            action = "add") {
  require(tidyverse)
  require(ComplexHeatmap)
  require(circlize)
  require(cluster)  # For silhouette scores
  require(vegan)  # For Gower's distance
  require(FactoMineR)  # Handling categorical data
  require(factoextra)
  require(dynamicTreeCut)
  require(densityClust)
  
  message("Starting clustering analysis...")
  
  # Validate exposure_cols input
  if (is.null(exposure_cols) || length(exposure_cols) == 0) {
    stop("Please specify `exposure_cols` as a character vector of `colData` column names.")
  }
  
  # Ensure selected exposure columns exist in colData
  available_variables <- intersect(exposure_cols, colnames(colData(expomicset)))
  if (length(available_variables) == 0) {
    stop("None of the specified `exposure_cols` are available in `colData`. Clustering cannot proceed.")
  }
  
  # Extract numeric exposure data for clustering
  exposure_data <- colData(expomicset) |> 
    as.data.frame() |> 
    dplyr::select(all_of(available_variables))
  
  # Remove missing/infinite values
  exposure_data <- exposure_data |> 
    na.omit() |> 
    mutate(across(everything(), ~ ifelse(is.infinite(.), NA, .))) |> 
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
    daisy(exposure_data, metric = "gower")
  } else {
    dist(exposure_data, method = dist_method)
  }
  
  # Function to determine optimal k based on the selected clustering approach
  determine_k <- function(dist_matrix, cluster_method) {
    if (clustering_approach == "diana") {
      sample_cluster <- diana(as.dist(dist_matrix))
      height_diffs <- diff(sample_cluster$height)
      cutoff_index <- which.max(height_diffs)
      return(length(sample_cluster$height) - cutoff_index)
      
    } else if (clustering_approach == "gap") {
      gap_stat <- clusGap(as.matrix(dist_matrix), FUN = hcut, K.max = 20, B = 50)
      return(maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"]))
      
    } else if (clustering_approach == "elbow") {
      fviz_nbclust(as.matrix(dist_matrix), FUN = hcut, method = "wss")
      return(3)  # Adjust manually if needed
      
    } else if (clustering_approach == "dynamic") {
      sample_cluster <- hclust(as.dist(dist_matrix), method = cluster_method)
      cut_clusters <- cutreeDynamic(dendro = as.dendrogram(sample_cluster), distM = as.matrix(dist_matrix), deepSplit = 2)
      return(length(unique(cut_clusters)))
      
    } else if (clustering_approach == "density") {
      dclust <- densityClust(dist_matrix, gaussian = TRUE)
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
   heatmap <- Heatmap(
    matrix = t(scale(exposure_data)),  # Scale and transpose for clustering
    name = "Scaled Exposure",
    cluster_columns = sample_cluster,
    top_annotation = HeatmapAnnotation(
      df=data.frame(
      Cluster = as.character(sample_groups),
      id = names(sample_groups)) |> 
        column_to_rownames("id")),
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8, angle = 45),
    heatmap_legend_param = list(title = "Exposure\nValue"),
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
  
  if (action == "add") {
    # Save clustering results in metadata
    metadata(expomicset)$sample_clustering <- list(
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
  
  # ta=HeatmapAnnotation(df=expom@metadata$sample_clustering$a |> 
  #                        as.data.frame() |> 
  #                        t() |> 
  #                        as.data.frame() |> 
  #                        rownames_to_column("id_to_map") |>
  #                        left_join(colData(expom) |> 
  #                                    as.data.frame() |> 
  #                                    dplyr::select(all_of(c("age","race"))) |> 
  #                                    rownames_to_column("id_to_map"),
  #                                  by = "id_to_map") |> 
  #                        column_to_rownames("id_to_map") |>
  #                        t() |> 
  #                        as.data.frame() |> 
  #                        (\(df) df[rownames(df) %in% c("age","race"),])() |> 
  #                        t() |> as.data.frame())
  # 
  # expom@metadata$sample_clustering$a |> 
  #   as.data.frame() |> 
  #   t() |> 
  #   as.data.frame() |> 
  #   rownames_to_column("id_to_map") |>
  #   left_join(colData(expom) |> 
  #               as.data.frame() |> 
  #               dplyr::select(all_of(c("age","race"))) |> 
  #               rownames_to_column("id_to_map"),
  #             by = "id_to_map") |> 
  #   column_to_rownames("id_to_map") |>
  #   t() |> 
  #   as.data.frame() |> 
  #   (\(df) df[!rownames(df) %in% c("age","race"),])() |> 
  #   mutate_all(as.numeric) |> 
  #   as.matrix() |> 
  #   Heatmap(  
  #     name = "Scaled Exposure",
  #     cluster_columns = expom@metadata$sample_clustering$sample_cluster,
  #     show_row_dend = TRUE,
  #     show_column_dend = TRUE,
  #     row_names_gp = gpar(fontsize = 8),
  #     column_names_gp = gpar(fontsize = 8, angle = 45),
  #     heatmap_legend_param = list(title = "Exposure\nValue"),
  #     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  #     top_annotation = ta
  #   )
}
