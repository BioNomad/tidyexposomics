#' Cluster Exposures Based on Scaled Exposure Data
#'
#' Performs hierarchical clustering of exposure variables using scaled
#' exposure data from `colData(exposomicset)` and records each exposure's
#' cluster in the codebook (`metadata(exposomicset)$codebook`).
#'
#' @param exposomicset A `MultiAssayExperiment` object containing omics
#' and exposure data.
#' @param exposure_cols A character vector of column names in
#' `colData(exposomicset)` to cluster.
#' @param dist_method A character string specifying the distance metric.
#' `"correlation"` (default) uses `1 - abs(cor)`, so exposures that are
#' strongly associated (positively or negatively) cluster together. Any
#' method accepted by `stats::dist()` can also be used.
#' @param user_k An integer specifying the number of clusters.
#' If `NULL`, an optimal `k` is determined.
#' @param cluster_method A character string specifying the hierarchical
#' clustering method. Default is `"ward.D"`.
#' @param clustering_approach A character string specifying the method
#' for determining `k` (`"diana"`, `"gap"`, `"elbow"`, `"dynamic"`,
#' or `"density"`). Default is `"diana"`.
#' @param codebook_var_col A character string naming the codebook column
#' that holds exposure names. Default is `"variable"`.
#' @param action A character string specifying `"add"`
#' (store results in metadata) or `"get"` (return clustering results).
#' Default is `"add"`.
#'
#' @return If `action="add"`, returns the updated `exposomicset` with an
#' `exposure_cluster` column in `metadata(exposomicset)$codebook`.
#' If `action="get"`, returns a list with:
#' \item{exposure_cluster}{A hierarchical clustering object (`hclust`).}
#' \item{exposure_groups}{A named vector of exposure cluster assignments.}
#'
#' @examples
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' mae <- run_cluster_exposures(
#'     exposomicset = mae,
#'     exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
#'     clustering_approach = "diana"
#' )
#'
#' @importFrom dplyr select all_of any_of left_join
#' @importFrom cluster daisy diana clusGap maxSE
#' @export
run_cluster_exposures <- function(exposomicset,
                                  exposure_cols = NULL,
                                  dist_method = "correlation",
                                  user_k = NULL,
                                  cluster_method = "ward.D",
                                  clustering_approach = "diana",
                                  codebook_var_col = "variable",
                                  action = "add") {
    message("Starting exposure clustering analysis...")
    .check_suggested(pkg = "dynamicTreeCut")
    .check_suggested(pkg = "densityClust")
    .check_suggested(pkg = "tidyselect")

    if (is.null(exposure_cols) || length(exposure_cols) == 0) {
        stop("Please specify `exposure_cols` as column names in `colData`.")
    }

    available_variables <- intersect(
        exposure_cols,
        colnames(MultiAssayExperiment::colData(exposomicset))
    )

    if (length(available_variables) < 2) {
        stop("At least two `exposure_cols` must be available in `colData`.")
    }

    # Scaled numeric exposures (samples x exposures), zero-variance removed
    scaled_exposures <- MultiAssayExperiment::colData(exposomicset) |>
        as.data.frame() |>
        dplyr::select(dplyr::all_of(available_variables)) |>
        dplyr::select(tidyselect::where(is.numeric)) |>
        na.omit() |>
        dplyr::select(tidyselect::where(~ sd(.x) > 0)) |>
        scale()

    if (ncol(scaled_exposures) < 2) {
        stop("Fewer than two usable numeric exposures after filtering.")
    }

    # Distance between exposures (columns)
    exposure_dist <- if (dist_method == "correlation") {
        as.dist(1 - abs(cor(scaled_exposures)))
    } else {
        dist(t(scaled_exposures), method = dist_method)
    }

    k_exposures <- if (is.null(user_k)) {
        .determine_k(
            dist_matrix = exposure_dist,
            cluster_method = cluster_method,
            clustering_approach = clustering_approach
        )
    } else {
        user_k
    }

    exposure_cluster <- hclust(exposure_dist, method = cluster_method)
    exposure_groups <- cutree(exposure_cluster, k = k_exposures)

    message("Optimal number of clusters for exposures: ", k_exposures)

    if (action == "get") {
        return(list(
            exposure_cluster = exposure_cluster,
            exposure_groups = exposure_groups
        ))
    } else if (action != "add") {
        stop("Invalid action specified. Use 'add' or 'get'.")
    }

    # Create a codebook if one doesn't exist yet
    if (is.null(MultiAssayExperiment::metadata(exposomicset)$codebook)) {
        MultiAssayExperiment::metadata(exposomicset)$codebook <-
            data.frame(names(exposure_groups)) |>
            stats::setNames(codebook_var_col)
    }

    # Add cluster membership to the codebook (replace on re-run)
    MultiAssayExperiment::metadata(exposomicset)$codebook <-
        MultiAssayExperiment::metadata(exposomicset)$codebook |>
        as.data.frame() |>
        dplyr::select(-dplyr::any_of("exposure_cluster")) |>
        dplyr::left_join(
            data.frame(
                exposure_to_map = names(exposure_groups),
                exposure_cluster = paste0("G", exposure_groups)
            ),
            by = stats::setNames("exposure_to_map", codebook_var_col)
        )

    # Save clustering object in metadata
    MultiAssayExperiment::metadata(exposomicset)$quality_control$exposure_clustering <- list(
        exposure_cluster = exposure_cluster,
        exposure_groups = exposure_groups
    )

    # Record analysis step
    MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
        MultiAssayExperiment::metadata(exposomicset)$summary$steps,
        list(
            run_cluster_exposures = list(
                timestamp = Sys.time(),
                params = list(
                    exposure_cols = exposure_cols,
                    dist_method = dist_method,
                    user_k = user_k,
                    cluster_method = cluster_method,
                    clustering_approach = clustering_approach
                ),
                notes = paste0(
                    "Optimal number of clusters for exposures: ",
                    k_exposures
                )
            )
        )
    )

    exposomicset
}
