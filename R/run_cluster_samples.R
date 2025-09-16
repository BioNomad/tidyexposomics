#' Cluster Samples Based on Exposure Data
#'
#' Performs hierarchical clustering of samples using exposure data from
#' `colData(exposomicset)`.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing omics
#' and exposure data.
#' @param exposure_cols A character vector of column names in
#' `colData(exposomicset)` to use for clustering.
#' @param dist_method A character string specifying the distance metric
#' (`"euclidean"`, `"gower"`, etc.). If `NULL`, it is automatically determined.
#' @param user_k An integer specifying the number of clusters.
#' If `NULL`, an optimal `k` is determined.
#' @param cluster_method A character string specifying the hierarchical
#' clustering method. Default is `"ward.D"`.
#' @param clustering_approach A character string specifying the method
#'  for determining `k` (`"diana"`, `"gap"`, `"elbow"`, `"dynamic"`,
#'  or `"density"`). Default is `"diana"`.
#' @param action A character string specifying `"add"`
#' (store results in metadata) or `"get"` (return clustering results).
#'  Default is `"add"`.
#'
#' @details
#' This function:
#' - Extracts **numeric exposure data** from `colData(exposomicset)`.
#' - Computes a **distance matrix** (`"gower"` for mixed data,
#' `"euclidean"` for numeric).
#' - Determines the **optimal number of clusters (`k`)** using the
#' specified method.
#' - Performs **hierarchical clustering** (`hclust`) and assigns
#' samples to clusters.
#' - Generates a **heatmap** of scaled exposure values.
#' - Stores results in `metadata(exposomicset)$sample_clustering`
#' when `action="add"`.
#'
#' @return If `action="add"`, returns the updated `exposomicset`.
#' If `action="get"`, returns a list with:
#' \item{sample_cluster}{A hierarchical clustering object (`hclust`).}
#' \item{sample_groups}{A named vector of sample cluster assignments.}
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # determine sample clusters
#' mae <- run_cluster_samples(
#'     exposomicset = mae,
#'     exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
#'     clustering_approach = "diana"
#' )
#'
#' @importFrom dplyr select mutate left_join all_of everything across
#' @importFrom tibble column_to_rownames
#' @importFrom cluster daisy diana clusGap maxSE
#' @importFrom factoextra hcut fviz_nbclust
#' @export
run_cluster_samples <- function(exposomicset,
                                exposure_cols = NULL,
                                dist_method = NULL,
                                user_k = NULL,
                                cluster_method = "ward.D",
                                clustering_approach = "diana",
                                action = "add") {
    message("Starting clustering analysis...")
    .check_suggested(pkg = "dynamicTreeCut")
    .check_suggested(pkg = "densityClust")

    # Validate exposure_cols input
    if (is.null(exposure_cols) || length(exposure_cols) == 0) {
        stop("Please specify `exposure_cols` is a a column name in `colData`.")
    }

    # Ensure selected exposure columns exist in colData
    available_variables <- intersect(
        exposure_cols,
        colnames(MultiAssayExperiment::colData(exposomicset))
    )

    if (length(available_variables) == 0) {
        stop("None of the specified `exposure_cols` are available in `colData`.")
    }

    # Extract numeric exposure data for clustering
    exposure_data <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        dplyr::select(dplyr::all_of(available_variables))

    # Remove missing/infinite values
    exposure_data <- exposure_data |>
        na.omit() |>
        dplyr::mutate(
            dplyr::across(
                dplyr::everything(),
                ~ ifelse(is.infinite(.),
                    NA,
                    .
                )
            )
        ) |>
        na.omit()

    if (nrow(exposure_data) == 0) {
        stop("All selected exposure data contains missing or infinite values.")
    }

    # Determine appropriate distance metric
    if (is.null(dist_method)) {
        if (all(vapply(exposure_data, is.numeric, FUN.VALUE = logical(1)))) {
            dist_method <- "euclidean" # Continuous data
        } else {
            dist_method <- "gower" # Mixed data types
        }
    }

    # Compute distance matrix
    sample_dist <- if (dist_method == "gower") {
        cluster::daisy(exposure_data, metric = "gower")
    } else {
        dist(exposure_data, method = dist_method)
    }

    if (!is.null(user_k)) {
        k_samples <- user_k
    } else {
        # Determine optimal k clusters
        k_samples <- .determine_k(
            dist_matrix = sample_dist,
            cluster_method = cluster_method,
            clustering_approach = clustering_approach
        )
    }

    # Perform hierarchical clustering
    sample_cluster <- hclust(as.dist(sample_dist),
        method = cluster_method
    )

    # Cut dendrogram using optimal k
    sample_groups <- cutree(sample_cluster, k = k_samples)

    message("Optimal number of clusters for samples: ", k_samples)

    # add in sample groups for samples to the colData
    col_data <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        (\(df){
            df$id_to_map <- rownames(df)
            df
        })() |>
        dplyr::left_join(
            data.frame(
                id_to_map = names(sample_groups),
                sample_group = paste0("Group_", as.character(sample_groups))
            ),
            by = "id_to_map"
        ) |>
        tibble::column_to_rownames("id_to_map")

    MultiAssayExperiment::colData(exposomicset) <- col_data

    if (action == "add") {
        # Save clustering results in metadata
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$quality_control$sample_clustering <- list(
            sample_cluster = sample_cluster,
            sample_groups = sample_groups
        )

        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

        # Add analysis steps taken to metadata
        step_record <- list(
            run_cluster_samples = list(
                timestamp = Sys.time(),
                params = list(
                    exposure_cols = exposure_cols,
                    dist_method = dist_method,
                    user_k = user_k,
                    cluster_method = cluster_method,
                    clustering_approach = clustering_approach
                ),
                notes = paste0("Optimal number of clusters for samples: ", k_samples)
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else if (action == "get") {
        return(list(
            sample_cluster = sample_cluster,
            sample_groups = sample_groups
        ))
    } else {
        stop("Invalid action specified. Use 'add' or 'get'.")
    }
}

# --- Determine K Clusters Function -------------
#' Determine Optimal Number of Clusters (k)
#'
#' Internal function to determine the optimal number of clusters based on a
#' selected clustering approach. Supports multiple methods such as gap statistic,
#' elbow method, DIANA, dynamic tree cut, and density-based clustering.
#'
#' @param dist_matrix A distance matrix (or object coercible to `dist`) representing
#'   pairwise sample distances.
#' @param cluster_method The hierarchical clustering method to use (e.g., "ward.D2",
#'   "average") if applicable.
#' @param clustering_approach The clustering approach used to determine number
#' of clusters.
#'
#' @return An integer indicating the estimated optimal number of clusters.
#' @keywords internal
#' @importFrom cluster diana clusGap maxSE
#' @importFrom factoextra hcut fviz_nbclust
#' @noRd
.determine_k <- function(
    dist_matrix,
    cluster_method,
    clustering_approach) {
    if (clustering_approach == "diana") {
        # Determine optimal k using the height difference method
        sample_cluster <- cluster::diana(as.dist(dist_matrix))
        height_diffs <- diff(sample_cluster$height)
        cutoff_index <- which.max(height_diffs)
        return(length(sample_cluster$height) - cutoff_index)
    } else if (clustering_approach == "gap") {
        # Determine optimal k using the gap statistic
        gap_stat <- cluster::clusGap(as.matrix(dist_matrix),
            FUN = factoextra::hcut,
            K.max = 20,
            B = 50
        )
        return(cluster::maxSE(
            gap_stat$Tab[, "gap"],
            gap_stat$Tab[, "SE.sim"]
        ))
    } else if (clustering_approach == "elbow") {
        # Determine optimal k using the elbow method
        wss_plot <- factoextra::fviz_nbclust(as.matrix(dist_matrix),
            FUN = factoextra::hcut,
            method = "wss"
        )

        # Identify the first significant drop and ensure it's a number we can use
        k_optimal <- which.min(diff(diff(wss_plot$data$y))) + 1
        if (is.na(k_optimal) || k_optimal < 2) k_optimal <- 3
        return(k_optimal)
    } else if (clustering_approach == "dynamic") {
        .check_suggested(pkg = "dynamicTreeCut")
        # Determine optimal k using dynamic tree cut
        sample_cluster <- hclust(as.dist(dist_matrix), method = cluster_method)
        cut_clusters <- dynamicTreeCut::cutreeDynamic(
            dendro = sample_cluster,
            distM = as.matrix(as.dist(dist_matrix)),
            deepSplit = 2
        )
        return(length(unique(cut_clusters)))
    } else if (clustering_approach == "density") {
        .check_suggested(pkg = "densityClust")
        # Determine optimal k using density-based clustering
        dclust <- densityClust::densityClust(as.dist(dist_matrix),
            gaussian = TRUE
        )
        dclust <- densityClust::findClusters(dclust,
            rho = quantile(dclust$rho, 0.90),
            delta = quantile(dclust$delta, 0.90)
        )
        return(length(unique(dclust$clusters)))
    } else {
        stop("Invalid clustering approach selected.")
    }
}
