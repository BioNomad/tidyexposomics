#' Run Multi-Omics Integration
#'
#' Performs multi-omics integration using one of several available methods: MOFA, MCIA,
#' RGCCA, or DIABLO. This function takes a `MultiAssayExperiment` object with two or more
#' assays and computes shared latent factors across omics layers.
#'
#' @param expomicset A `MultiAssayExperiment` object with at least two assays.
#' @param method Character. Integration method to use. Options are `"MOFA"`, `"MCIA"`,
#'   `"RGCCA"`, or `"DIABLO"`.
#' @param n_factors Integer. Number of latent factors/components to compute. Default is 10.
#' @param scale Logical. Whether to scale each assay before integration. Default is `TRUE`.
#' @param outcome Character. Required if `method = "DIABLO"`. Name of outcome variable in `colData` used for supervised integration.
#' @param action Character. Whether to `"add"` results to the metadata or `"get"` them as a list. Default is `"add"`.
#'
#' @return If `action = "add"`, returns a `MultiAssayExperiment` with integration results
#' stored in `metadata(expomicset)$multiomics_integration$integration_results`. If `action = "get"`, returns a list with integration `method` and `result`.
#'
#' @details
#' - `"MOFA"` runs Multi-Omics Factor Analysis using the `MOFA2` package and returns a trained model.
#' - `"MCIA"` runs multi-co-inertia analysis using the `nipalsMCIA` package.
#' - `"RGCCA"` runs Regularized Generalized Canonical Correlation Analysis using the `RGCCA` package.
#' - `"DIABLO"` performs supervised integration using the `mixOmics` package and a specified outcome.
#'
#' @seealso \code{\link{plot_factor_scores}}, \code{\link{plot_top_factor_features}},
#'   \code{\link{run_factor_overlap}}, \code{\link{run_enrichment}}
#'
#' @examples
#' \dontrun{
#' expomicset <- run_multiomics_integration(
#'   expomicset,
#'   method = "MOFA",
#'   n_factors = 5,
#'   scale = TRUE
#' )
#'
#' expomicset <- run_multiomics_integration(
#'   expomicset,
#'   method = "DIABLO",
#'   n_factors = 3,
#'   outcome = "asthma_status"
#' )
#' }
#'
#' @export
run_multiomics_integration <- function(expomicset,
                                   method = "MCIA",
                                   n_factors = 10,
                                   scale=TRUE,
                                   outcome=NULL,
                                   action="add") {

  if(length(MultiAssayExperiment::experiments(expomicset)) < 2){
    stop("Multi-Omics Integration requires at least two assays in the MultiAssayExperiment object.")
  }

  if(scale){
    expomicset_mo <- .scale_multiassay(expomicset)
  }
  else{
    expomicset_mo <- expomicset
  }

  message("Running multi-omics integration using ", method, "...")

  result <- NULL

  # MOFA integration
  if (method == "MOFA") {
    message("Applying MOFA+ integration.")

    # Create MOFA object from the MultiAssayExperiment
    mofa <- MOFA2::create_mofa(expomicset_mo)

    # Set MOFA options
    model_opts <- MOFA2::get_default_model_options(mofa)
    model_opts$num_factors <- n_factors
    data_opts <- MOFA2::get_default_data_options(mofa)
    train_opts <- MOFA2::get_default_training_options(mofa)

    # Prepare & train MOFA model
    mofa <- MOFA2::prepare_mofa(
      object = mofa,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )

    outfile = file.path(tempdir(), "mofa_model.hdf5")
    mofa_trained <- MOFA2::run_mofa(mofa, outfile, use_basilisk = TRUE)

    # Load trained MOFA model
    result <- MOFA2::load_model(outfile)

    # NIPALS MCIA INTEGRATION
  } else if (method == "MCIA") {
    message("Applying MCIA with NIPALS.")

    # Run NIPALS MCIA on the MultiAssayExperiment

    result <- nipalsMCIA::nipals_multiblock(
      expomicset_mo,
      col_preproc_method = "colprofile",
      num_PCs = n_factors,
      tol = 1e-12,
      plots = "none")

  } else if (method == "RGCCA") {
    message("Applying RGCCA integration.")

    # Prepare data: list of assays (samples x features)
    x <- purrr::map(
      MultiAssayExperiment::experiments(expomicset_mo),
      ~ t(SummarizedExperiment::assay(.x))  # transpose to samples x features
    )

    # Use default shrinkage or manual choice (e.g., 1 for PCA-like)
    tau <- rep(1, length(x))  # no regularization

    # Run RGCCA
    rgcca_result <- RGCCA::rgcca(
      blocks = x,
      connection = 1 - diag(length(x)),
      tau = tau,
      ncomp = n_factors,
      scheme = "centroid",
      scale = FALSE,
      verbose = FALSE
    )

    result <- rgcca_result
  } else if (method == "DIABLO") {
    message("Applying DIABLO supervised integration.")

    if (is.null(outcome)) {
      stop("DIABLO requires an outcome variable for supervision.")
    }
    # Ensure that only the common samples are taken
    common <- purrr::map(
      MultiAssayExperiment::experiments(expomicset_mo),
      ~ colnames(.x)
    ) |>
      (\(lst) Reduce(intersect,lst))()

    # Subset the samples
    expomicset_mo <- expomicset_mo[,common]

    # Extract colData and outcome
    y <- MultiAssayExperiment::colData(expomicset_mo)[[outcome]]
    if (!is.factor(y)) y <- as.factor(y)

    # Prepare blocks (assays) as list of matrices (samples × features)
    blocks <- purrr::map(
      MultiAssayExperiment::experiments(expomicset_mo),
      ~ t(SummarizedExperiment::assay(.x))
    )

    # Use 1 component by default for each block
    ncomp <- rep(n_factors, length(blocks))
    names(ncomp) <- names(blocks)

    # Run DIABLO
    design <- matrix(1, ncol = length(blocks), nrow = length(blocks))
    diag(design) <- 0
    rownames(design) <- colnames(design) <- names(blocks)

    result <- mixOmics::block.splsda(
      X = blocks,
      Y = y,
      ncomp = n_factors,
      design = design
    )
  } else {
    stop("Invalid method. Choose from 'MOFA', 'MCIA', 'RGCCA', or 'DIABLO'.")
  }

  if(action=="add"){
    # Store results in MultiAssayExperiment metadata
    MultiAssayExperiment::metadata(expomicset)$multiomics_integration$integration_results <- list(
      method = method,
      result = result
    )

    # Add analysis step record
    step_record <- list(
      run_multiomics_integration = list(
        timestamp = Sys.time(),
        params = list(
          method = method,
          n_factors = n_factors,
          scale = scale
        ),
        notes = paste0(
          "Performed multi-omics integration using ", method,
          " with ", n_factors, " latent factors. Scaling was ",
          ifelse(scale, "enabled.", "disabled.")
        )
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )


    return(expomicset)
  }else if (action=="get"){
    return(list(
      method = method,
      result = result
    ))
  }else{
    stop("Invalid action. Choose from 'add' or 'get'.")
  }
}


# else if (method == "MCCA") {
#   # MCCA integration using PMA package
#   x <- purrr::map(
#     .x = MultiAssayExperiment::experiments(expomicset_mo),
#     .f = ~ SummarizedExperiment::assay(.x)
#   )
#   # Ensure all assays have the same columns
#   x <- purrr::map(
#     .x = x,
#     .f = ~ t(.x[, Reduce(intersect,lapply(x, colnames))])
#   )
#
#   # Run MultiCCA
#   y <- PMA::MultiCCA(
#     xlist = x,
#     ncomponents = n_factors)
#
#   # Extract results
#   y_ws <- y$ws
#
#   # Match the names of the results to the original assays
#   names(y_ws) <- names(x)
#
#   # Match rownames
#   y_ws <- purrr::imap(
#     .x = y_ws,
#     .f = ~ {
#       rownames(.x) <- colnames(x[[.y]])
#       .x
#     }
#   )
#
#   # Grab sample level scores
#   sample_scores <- purrr::map2(
#     .x = x,         # data: samples x features
#     .y = y$ws,      # weights: features x components
#     .f = ~ .x %*% .y   # scores: samples x components
#   )
#
#   # return result
#   result <- list(
#     weights = y_ws,
#     sample_scores = sample_scores
#   )
#
# }
