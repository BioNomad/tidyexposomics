#' Calculate Exposure Impact Based on DEG Network Centrality
#'
#' This function computes a network-based measure of exposure impact by integrating
#' differential abundance (DA) results, DA-exposure correlations, and network centrality.
#' It returns either a summary of exposure impact statistics or stores them in the
#' `MultiAssayExperiment` object's metadata.
#'
#' @param expomicset A \code{MultiAssayExperiment} object containing prior results
#'   from sensitivity analysis and exposure correlation analysis.
#' @param robust Logical. If \code{TRUE}, checks that sensitivity analysis has been run
#'   and uses robust (stability-based) features. Default is \code{TRUE}.
#' @param stability_threshold Numeric threshold used to filter features by stability
#'   score. If \code{NULL}, the threshold stored in metadata is used.
#' @param pval_col Character. Column name in the DA results indicating the adjusted
#'   p-value. Default is \code{"adj.P.Val"}.
#' @param pval_thresh Numeric. Threshold used to define differentially abundant features.
#'   Default is \code{0.1}.
#' @param action Character. If \code{"add"}, stores results in metadata; if \code{"get"},
#'   returns results as a list. Must be either \code{"add"} or \code{"get"}.
#'
#' @return Either:
#'   \itemize{
#'     \item The original \code{MultiAssayExperiment} object with added exposure impact metadata (if \code{action = "add"}), or
#'     \item A named list with:
#'       \describe{
#'         \item{exposure_impact_degree}{Data frame of feature centrality statistics by exposure.}
#'         \item{exposure_impact_deg_n}{Summary of how many DEGs are associated with each exposure.}
#'       }
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Filters stable features (using stability scores if \code{robust = TRUE}).
#'   \item Runs exposure correlation analysis on those features.
#'   \item Builds a network of correlated features and calculates node centrality.
#'   \item Merges centrality with DEG-exposure associations.
#'   \item Summarizes exposure impact using average centrality and DEG counts.
#' }
#'
#' @seealso \code{\link{run_sensitivity_analysis}}, \code{\link{correlate_exposoures_degs}}, \code{\link{correlate_exposoures_omics}}
#'
#'
#' @export

run_exposure_impact <- function(
    expomicset,
    robust=TRUE,
    stability_threshold=NULL,
    pval_col="adj.P.Val",
    pval_thresh=0.1,
    action="add"
){
  # This function will take a multiassay experiment object
  # search the meta data
  # identify differential abundance results and da-exposure correlations
  # create a network based on correlation of the differential abundance results
  # and then get the degree of those correlations
  # then map back to the da-exposure correlations
  # and identify the degree for features that
  # are associated with a particular exposure
  # and then plot the degree of those features

  require(ggplot2)

  # Check to see if the stability analysis has been run
  if(robust & !("sensitivity_analysis" %in% names(MultiAssayExperiment::metadata(expomicset)))) {
    stop("Please run `run_senstivity_analysis()` first.")
  }

  # Check to see if the deg correlation analysis has been run
  if(!("omics_exposure_deg_correlation" %in% names(MultiAssayExperiment::metadata(expomicset)))) {
    stop("Please run `correlate_exposoures_degs()` first.")
  }

  if(!is.null(stability_threshold)) {
    stability_threshold <- stability_threshold
  } else{
    stability_threshold <- expomicset |>
      MultiAssayExperiment::metadata() |>
      purrr::pluck("sensitivity_analysis") |>
      purrr::pluck("score_thresh")
  }

  da_res <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("sensitivity_analysis") |>
    purrr:::pluck("feature_stability") |>
    dplyr::filter(stability_score>stability_threshold) |>
    dplyr::select(feature,exp_name) |>
    dplyr::distinct()

  expomicset_tmp <- expomicset |>
    correlate_exposoures_omics(
      variable_map = da_res |>
        dplyr::ungroup() |>
        dplyr::rename(variable=feature),
      action = "add"
    )

  node_centrality <- expomicset_tmp |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("omics_exposure_correlation") |>
    tidygraph::as_tbl_graph() |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(centrality=tidygraph::centrality_degree()) |>
    tidygraph::as_tibble() |>
    as.data.frame() |>
    dplyr::rename(feature=name)

  exposure_impact_degree <- expomicset_tmp |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("omics_exposure_deg_correlation") |>
    dplyr::left_join(node_centrality,by = "feature") |>
    dplyr::filter(!is.na(centrality)) |>
    dplyr::group_by(exposure) |>
    dplyr::mutate(
      mean=mean(centrality,na.rm = TRUE),
      n_exp=dplyr::n()
    ) |>
    dplyr::ungroup()

  exposure_impact_deg_n <- expomicset_tmp |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("differential_abundance") |>
    dplyr::filter(!!dplyr::sym(pval_col)<pval_thresh) |>
    dplyr::left_join(expomicset_tmp |>
                       MultiAssayExperiment::metadata() |>
                       purrr::pluck("omics_exposure_deg_correlation") |>
                dplyr::select(feature,exp_name) |>
                dplyr::distinct() |>
                dplyr::mutate(exposure_association = "yes"),
              by = c("feature", "exp_name")
    ) |>
    dplyr::mutate(exposure_association=ifelse(is.na(exposure_association),
                                       "no",
                                       exposure_association)) |>
    dplyr::group_by(exp_name) |>
    dplyr::summarise(
      total_deg=dplyr::n(),
      n_assoc_exp=length(exposure_association[exposure_association == "yes"]),
      percent_assoc_exp = n_assoc_exp/total_deg*100
    )

  exposure_impact <- list(
    exposure_impact_degree=exposure_impact_degree,
    exposure_impact_deg_n=exposure_impact_deg_n
  )

  if(action=="add"){
    # Store selected features
    MultiAssayExperiment::metadata(expomicset)$exposure_impact <- exposure_impact

    # Add analysis steps taken to metadata
    MultiAssayExperiment::metadata(expomicset)$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$steps,
      "run_exposure_impact"
    )

    return(expomicset)
  }else if (action=="get"){
    return(exposure_impact)
  }else{
    stop("Invalid action. Choose 'add' or 'get'.")
  }


}
