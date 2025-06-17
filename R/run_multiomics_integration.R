#' Perform Multi-Omics Integration
#'
#' Integrates multiple assays in a `MultiAssayExperiment` object using either MOFA+ or MCIA.
#'
#' @param expomicset A `MultiAssayExperiment` object containing at least two omics assays.
#' @param method A character string specifying the integration method. Options: `"MOFA"`, `"MCIA"`. Default is `"MCIA"`.
#' @param n_factors An integer specifying the number of latent factors/components to compute. Default is `10`.
#' @param scale A logical indicating whether to standardize (`Z-score`) the assays before integration. Default is `TRUE`.
#' @param action A character string specifying whether to store (`"add"`) or return (`"get"`) the results. Default is `"add"`.
#'
#' @details
#' This function:
#' - **MOFA+ (`method = "MOFA"`)**:
#'   - Creates a MOFA model from `expomicset`.
#'   - Trains the model using default MOFA parameters.
#'   - Saves and loads the trained model.
#' - **MCIA (`method = "MCIA"`)**:
#'   - Applies **NIPALS MCIA** via `nipalsMCIA::nipals_multiblock()`.
#' - **Scaling**: If `scale = TRUE`, all assays are standardized using `.scale_multiassay()`.
#' - **Output Handling**:
#'   - `"add"`: Stores results in `metadata(expomicset)$integration_results`.
#'   - `"get"`: Returns a list containing `method` and `result` (trained model).
#'
#' @return A `MultiAssayExperiment` object with integration results added to metadata (if `action = "add"`) or a list with:
#' \item{method}{The integration method used (`"MOFA"` or `"MCIA"`).}
#' \item{result}{The trained integration model.}
#'
#' @examples
#' \dontrun{
#' expom <- run_multiomics_integration(
#'   expomicset = expom,
#'   method = "MCIA",
#'   n_factors = 10,
#'   scale = TRUE
#' )
#' }
#'
#' @export
run_multiomics_integration <- function(expomicset,
                                   method = "MCIA",
                                   n_factors = 10,
                                   scale=TRUE,
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
    message("Applying MOFA+ integration...")

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
    message("Applying MCIA with NIPALS...")

    # Run NIPALS MCIA on the MultiAssayExperiment

    result <- nipalsMCIA::nipals_multiblock(
      expomicset_mo,
      col_preproc_method = "colprofile",
      num_PCs = n_factors,
      tol = 1e-12,
      plots = "none")

  } else if (method == "MCCA") {
    # MCCA integration using PMA package
    x <- purrr::map(
      .x = MultiAssayExperiment::experiments(expomicset_mo),
      .f = ~ SummarizedExperiment::assay(.x)
    )
    # Ensure all assays have the same columns
    x <- purrr::map(
      .x = x,
      .f = ~ t(.x[, Reduce(intersect,lapply(x, colnames))])
    )

    # Run MultiCCA
    y <- PMA::MultiCCA(
      xlist = x,
      ncomponents = n_factors)

    # Extract results
    y_ws <- y$ws

    # Match the names of the results to the original assays
    names(y_ws) <- names(x)

    # Match rownames
    y_ws <- purrr::imap(
      .x = y_ws,
      .f = ~ {
        rownames(.x) <- colnames(x[[.y]])
        .x
      }
    )

    # Grab sample level scores
    sample_scores <- purrr::map2(
      .x = x,         # data: samples x features
      .y = y$ws,      # weights: features x components
      .f = ~ .x %*% .y   # scores: samples x components
    )

    # return result
    result <- list(
      weights = y_ws,
      sample_scores = sample_scores
    )

  } else{
    stop("Invalid method. Choose from 'MOFA', 'MCIA', or 'MCCA'.")
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
