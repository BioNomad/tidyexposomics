#' Cluster Samples Based on Exposure Data
#'
#' Performs hierarchical clustering of samples using exposure data from `colData(expomicset)`.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics and exposure data.
#' @param exposure_cols A character vector of column names in `colData(expomicset)` to use for clustering.
#' @param dist_method A character string specifying the distance metric (`"euclidean"`, `"gower"`, etc.). If `NULL`, it is automatically determined.
#' @param user_k An integer specifying the number of clusters. If `NULL`, an optimal `k` is determined.
#' @param cluster_method A character string specifying the hierarchical clustering method. Default is `"ward.D"`.
#' @param clustering_approach A character string specifying the method for determining `k` (`"diana"`, `"gap"`, `"elbow"`, `"dynamic"`, or `"density"`). Default is `"diana"`.
#' @param action A character string specifying `"add"` (store results in metadata) or `"get"` (return clustering results). Default is `"add"`.
#'
#' @details
#' This function:
#' - Extracts **numeric exposure data** from `colData(expomicset)`.
#' - Computes a **distance matrix** (`"gower"` for mixed data, `"euclidean"` for numeric).
#' - Determines the **optimal number of clusters (`k`)** using the specified method.
#' - Performs **hierarchical clustering** (`hclust`) and assigns samples to clusters.
#' - Generates a **heatmap** of scaled exposure values.
#' - Stores results in `metadata(expomicset)$sample_clustering` when `action="add"`.
#'
#' @return If `action="add"`, returns the updated `expomicset`.
#' If `action="get"`, returns a list with:
#' \item{sample_cluster}{A hierarchical clustering object (`hclust`).}
#' \item{sample_groups}{A named vector of sample cluster assignments.}
#' \item{heatmap}{A `ComplexHeatmap` object visualizing sample clustering.}
#'
#' @examples
#' \dontrun{
#' expom <- run_cluster_samples(
#'   expomicset = expom,
#'   exposure_cols = c("PM2.5", "NO2"),
#'   dist_method = "gower",
#'   clustering_approach = "gap"
#' )
#' }
#'
#' @export
run_cluster_samples <- function(expomicset,
                            exposure_cols = NULL,
                            dist_method = NULL,
                            user_k = NULL,
                            cluster_method = "ward.D",
                            clustering_approach = "diana",
                            action = "add") {

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
      gap_stat <- cluster::clusGap(as.matrix(dist_matrix), FUN = factoextra::hcut, K.max = 20, B = 50)
      return(cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"]))

    } else if (clustering_approach == "elbow") {

      # Determine optimal k using the elbow method
      wss_plot <- factoextra::fviz_nbclust(as.matrix(dist_matrix), FUN = factoextra::hcut, method = "wss")

      # Identify the first significant drop and ensure it's a number we can use
      k_optimal <- which.min(diff(diff(wss_plot$data$y))) + 1
      if (is.na(k_optimal) || k_optimal < 2) k_optimal <- 3
      return(k_optimal)

    } else if (clustering_approach == "dynamic") {

      # Determine optimal k using dynamic tree cut
      sample_cluster <- hclust(as.dist(dist_matrix), method = cluster_method)
      cut_clusters <- dynamicTreeCut::cutreeDynamic(
        dendro = sample_cluster,
        distM = as.matrix(as.dist(dist_matrix)),
        deepSplit = 2)
      return(length(unique(cut_clusters)))

    } else if (clustering_approach == "density") {
      # Determine optimal k using density-based clustering
      dclust <- densityClust::densityClust(as.dist(dist_matrix), gaussian = TRUE)
      dclust <- densityClust::findClusters(dclust, rho = quantile(dclust$rho, 0.90), delta = quantile(dclust$delta, 0.90))
      return(length(unique(dclust$clusters)))

    } else {
      stop("Invalid clustering approach selected.")
    }
  }

  if(!is.null(user_k)) {
    k_samples <- user_k
  } else {
    # Determine optimal k clusters
    k_samples <- determine_k(sample_dist, cluster_method)
  }

  # Perform hierarchical clustering
  sample_cluster <- hclust(as.dist(sample_dist), method = cluster_method)

  # Cut dendrogram using optimal k
  sample_groups <- cutree(sample_cluster, k = k_samples)

  message("Optimal number of clusters for samples: ", k_samples)

  # add in sample groups for samples to the colData
  col_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    (\(df){df$id_to_map=rownames(df);df})() |>
    dplyr::left_join(data.frame(
      id_to_map = names(sample_groups),
      sample_group = paste0("Group_",as.character(sample_groups))
    ),
                     by="id_to_map") |>
    tibble::column_to_rownames("id_to_map")

  MultiAssayExperiment::colData(expomicset) <- col_data

  # # Generate heatmap
  #  heatmap <- ComplexHeatmap::Heatmap(
  #   matrix = t(scale(exposure_data)),
  #   name = "Scaled Exposure",
  #   cluster_columns = sample_cluster,
  #   top_annotation = ComplexHeatmap::HeatmapAnnotation(
  #     df=data.frame(
  #     Cluster = as.character(sample_groups),
  #     id = names(sample_groups)) |>
  #       tibble::column_to_rownames("id")),
  #   show_row_dend = TRUE,
  #   show_column_dend = TRUE,
  #   row_names_gp = grid::gpar(fontsize = 8),
  #   column_names_gp = grid::gpar(fontsize = 8, angle = 45),
  #   heatmap_legend_param = list(title = "Exposure\nValue"),
  #   col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # )

  if (action == "add") {
    # Save clustering results in metadata
    MultiAssayExperiment::metadata(expomicset)$sample_clustering <- list(
      sample_cluster = sample_cluster,
      sample_groups = sample_groups
    )
    return(expomicset)
  } else if (action == "get") {
    return(list(
      sample_cluster = sample_cluster,
      sample_groups = sample_groups
    ))
  } else {
    stop("Invalid action specified. Use 'add' or 'get'.")
  }
}
