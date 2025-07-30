#' Filter low-quality features in omics assays
#'
#' This function applies variance- or expression-based filtering
#' across one or more assays within a `MultiAssayExperiment` object.
#' It is useful for removing low-quality or uninformative features
#'  before downstream analysis.
#'
#' @param expomicset A `MultiAssayExperiment` object containing omics assays.
#' @param method Filtering method: either `"variance"` or `"expression"`.
#' @param assays Character vector of assay names to filter.
#' If `NULL`, all assays are filtered.
#' @param assay_name Name or index of the assay within each
#' `SummarizedExperiment` to use.
#' @param min_var Minimum variance threshold (used if `method = "variance"`).
#' @param min_value Minimum expression value (used if `method = "expression"`).
#' @param min_prop Minimum proportion of samples exceeding `min_value`
#' (used if `method = "expression"`).
#' @param verbose Whether to print messages for each assay being filtered.
#'
#' @return A filtered `MultiAssayExperiment` object with
#' updated assays and step record.
#' @examples
#' # Filter the proteomics assay by variance
#' filtered_mae <-  filter_omics(
#'   expomicset = make_example_data(return_mae=TRUE),
#'   method = c("variance"),
#'   assays = "proteomics",
#'   assay_name = 1,
#'   min_var = 0.01,
#'   verbose = TRUE)
#'
#'
#' @export
filter_omics <- function(expomicset,
                         method = c("variance", "expression"),
                         assays = NULL,
                         assay_name = 1,
                         min_var = 1e-5,
                         min_value = 5,
                         min_prop = 0.7,
                         verbose = TRUE) {

  method <- match.arg(method)
  if (!inherits(expomicset, "MultiAssayExperiment")) {
    stop("Input must be a MultiAssayExperiment object.")
  }

  exp_list <- MultiAssayExperiment::experiments(expomicset)
  assay_names <- names(exp_list)
  if (is.null(assays)) assays <- assay_names
  assays <- intersect(assays, assay_names)

  removed_features <- list()

  # Now filter_single RETURNS a list: list(se = ..., stats = ...)
  filter_single <- function(se, name) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    total <- nrow(mat)

    keep <- switch(method,
                   variance = {
                     row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
                     row_vars > min_var
                   },
                   expression = {
                     prop_expressed <- rowMeans(mat > min_value, na.rm = TRUE)
                     prop_expressed >= min_prop
                   }
    )

    removed <- sum(!keep)
    stats <- list(
      removed = removed,
      total = total,
      kept = total - removed
    )

    if (verbose) {
      message(sprintf("Filtered %d of %d features from '%s' using method '%s'",
                      removed, total, name, method))
    }

    list(se = se[keep, , drop = FALSE], stats = stats)
  }

  filtered_exps <- list()
  for (name in assay_names) {
    se <- exp_list[[name]]
    if (name %in% assays) {
      if (verbose) message(sprintf("Filtering assay: %s", name))
      result <- filter_single(se, name)
      filtered_exps[[name]] <- result$se
      removed_features[[name]] <- result$stats
    } else {
      filtered_exps[[name]] <- se
    }
  }

  filtered_exps_list <- MultiAssayExperiment::ExperimentList(filtered_exps)
  MultiAssayExperiment::experiments(expomicset) <- filtered_exps_list

  # Update metadata with step records
  for (assay in names(removed_features)) {
    info <- removed_features[[assay]]
    step_name <- paste0("filter_omics_", assay)
    step_record <- list()
    step_record[[step_name]] <- list(
      timestamp = Sys.time(),
      params = list(
        method = method,
        assay = assay,
        assay_name = assay_name,
        min_var = min_var,
        min_value = min_value,
        min_prop = min_prop
      ),
      notes = sprintf(
        "Filtered omics features from '%s'
        Using method = '%s': %d removed of %d (%.1f%%).",
        assay,
        method,
        info$removed,
        info$total,
        100 * info$removed / info$total
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )
  }

  return(expomicset)
}

# filter_omics <- function(expomicset,
#                          method = c("variance", "expression"),
#                          assays = NULL,
#                          assay_name = 1,
#                          min_var = 1e-5,
#                          min_value = 5,
#                          min_prop = 0.7,
#                          verbose = TRUE) {
#   method <- match.arg(method)
#
#   if (!inherits(expomicset, "MultiAssayExperiment")) {
#     stop("Input must be a MultiAssayExperiment object.")
#   }
#
#   exp_list <- MultiAssayExperiment::experiments(expomicset)
#   assay_names <- names(exp_list)
#
#   if (is.null(assays)) assays <- assay_names
#   assays <- intersect(assays, assay_names)
#
#   removed_features <- list()
#
#   filter_single <- function(se, name) {
#     mat <- SummarizedExperiment::assay(se, assay_name)
#     total <- nrow(mat)
#
#     keep <- switch(method,
#                    variance = {
#                      row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
#                      row_vars > min_var
#                    },
#                    expression = {
#                      prop_expressed <- rowMeans(mat > min_value, na.rm = TRUE)
#                      prop_expressed >= min_prop
#                    })
#
#     removed <- sum(!keep)
#     removed_features[[name]] <<- list(
#       removed = removed,
#       total = total,
#       kept = total - removed
#     )
#
#     if (verbose) {
#       message(sprintf("Filtered %d of %d features from '%s' using method '%s'",
#                       removed, total, name, method))
#     }
#
#     se[keep, , drop = FALSE]
#   }
#
#   filtered_exps <- lapply(assay_names, function(name) {
#     se <- exp_list[[name]]
#     if (name %in% assays) {
#       if (verbose) message(sprintf("Filtering assay: %s", name))
#       filter_single(se, name)
#     } else {
#       se
#     }
#   })
#
#   names(filtered_exps) <- assay_names
#   MultiAssayExperiment::experiments(expomicset) <- MultiAssayExperiment::ExperimentList(filtered_exps)
#
#   # Add step records per assay
#   for (assay in names(removed_features)) {
#     info <- removed_features[[assay]]
#     step_name <- paste0("filter_omics_", assay)
#     step_record <- list()
#     step_record[[step_name]] <- list(
#       timestamp = Sys.time(),
#       params = list(
#         method = method,
#         assay = assay,
#         assay_name = assay_name,
#         min_var = min_var,
#         min_value = min_value,
#         min_prop = min_prop
#       ),
#       notes = sprintf(
#         "Filtered omics features from '%s' using method = '%s': %d removed of %d (%.1f%%).",
#         assay, method, info$removed, info$total, 100 * info$removed / info$total
#       )
#     )
#
#     MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
#       MultiAssayExperiment::metadata(expomicset)$summary$steps,
#       step_record
#     )
#   }
#
#   return(expomicset)
# }
#
