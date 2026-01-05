#' Run Multi-Omics Integration
#'
#' Performs multi-omics integration using one of several available methods:
#'  MOFA, MCIA, RGCCA, or DIABLO. This function takes a `MultiAssayExperiment`
#'  object with two or more assays and computes shared latent factors
#'  across omics layers.
#'
#' @param exposomicset A `MultiAssayExperiment` object with at least two assays.
#' @param method Character. Integration method to use. Options are
#' `"MOFA"`, `"MCIA"`, `"RGCCA"`, or `"DIABLO"`.
#' @param n_factors Integer. Number of latent factors/components to compute.
#' Default is 10.
#' @param scale Logical. Whether to scale each assay before integration.
#'  Default is `TRUE`.
#' @param outcome Character. Required if `method = "DIABLO"`.
#' Name of outcome variable in `colData` used for supervised integration.
#' @param max.iter numeric. Option to increase the number of iterations for
#'  `mixOmics::block.splsda` the default is 500.
#' @param near.zero.var Logical. Option to remove variables with near zero
#'  variance for `mixOmics::block.splsda`, default is `TRUE` .
#' @param action Character. Whether to `"add"` results to the metadata or
#' `"get"` them as a list. Default is `"add"`.
#'
#' @return If `action = "add"`, returns a `MultiAssayExperiment` with
#'  integration results
#' stored in `metadata(exposomicset)$multiomics_integration$integration_results`.
#'  If `action = "get"`, returns a list with integration `method` and `result`.
#'
#' @details
#' - `"MOFA"` runs Multi-Omics Factor Analysis using the `MOFA2` package and
#' returns a trained model.
#' - `"MCIA"` runs multi-co-inertia analysis using the `nipalsMCIA` package.
#' - `"RGCCA"` runs Regularized Generalized Canonical Correlation Analysis
#'  using the `RGCCA` package.
#' - `"DIABLO"` performs supervised integration using the `mixOmics` package
#' and a specified outcome.
#'
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # perform multiomics integration
#' mae <- run_multiomics_integration(
#'     mae,
#'     method = "DIABLO",
#'     outcome = "smoker",
#'     n_factors = 3
#' )
#'
#' @export
run_multiomics_integration <- function(
  exposomicset,
  method = "MCIA",
  n_factors = 10,
  scale = TRUE,
  outcome = NULL,
  max.iter = 500,
  near.zero.var = TRUE,
  action = "add"
) {
    if (length(MultiAssayExperiment::experiments(exposomicset)) < 2) {
        stop("Multi-Omics Integration requires at least two assays.")
    }

    if (scale) {
        exposomicset_mo <- .scale_multiassay(exposomicset)
    } else {
        exposomicset_mo <- exposomicset
    }

    message("Running multi-omics integration using ", method, "...")

    result <- switch(method,
        "MOFA" = .run_mofa2(
            exposomicset_mo = exposomicset_mo,
            n_factors = n_factors
        ),
        "MCIA" = .run_mcia(
            exposomicset_mo = exposomicset_mo,
            n_factors = n_factors
        ),
        "RGCCA" = .run_rgcca(
            exposomicset_mo = exposomicset_mo,
            n_factors = n_factors
        ),
        "DIABLO" = .run_diablo(
            exposomicset_mo = exposomicset_mo,
            n_factors = n_factors,
            outcome = outcome,
            max.iter = max.iter,
            near.zero.var = near.zero.var
        )
    )

    if (action == "add") {
        # Store results in MultiAssayExperiment metadata
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$multiomics_integration$integration_results <- list(
            method = method,
            result = result
        )
        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

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

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )


        return(exposomicset)
    } else if (action == "get") {
        return(list(
            method = method,
            result = result
        ))
    } else {
        stop("Invalid action. Choose from 'add' or 'get'.")
    }
}

# --- MOFA2 Integration Function ---------------
#' Decide which backend to use (basilisk vs reticulate)
#'
#' @return "basilisk" or "reticulate"
#' @keywords internal
#' @noRd
.mofa_backend <- function() {
    # user override wins
    backend <- getOption("tidyexposomics.mofa.backend", NULL)
    if (!is.null(backend)) {
        return(backend)
    }

    # auto-detect: Apple Silicon, reticulate, otherwise basilisk
    if (Sys.info()[["sysname"]] == "Darwin" &&
        grepl("arm64", R.version$platform)) {
        return("reticulate")
    } else {
        return("basilisk")
    }
}

#' Check that the chosen backend is safe to use
#' @keywords internal
#' @noRd
.check_mofa_safe <- function(backend = .mofa_backend()) {
    if (identical(backend, "reticulate")) {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("reticulate backend requested, but 'reticulate' is not installed.", call. = FALSE)
        }
        ret <- asNamespace("reticulate")
        if (!isTRUE(ret$py_available(initialize = FALSE))) {
            stop("reticulate backend requested, but no Python is available.", call. = FALSE)
        }
        if (!isTRUE(ret$py_module_available("mofapy2"))) {
            stop("reticulate backend requested, but Python module 'mofapy2' is missing.", call. = FALSE)
        }
    } else if (identical(backend, "basilisk")) {
        if (!requireNamespace("MOFA2", quietly = TRUE)) {
            stop("basilisk backend requested, but 'MOFA2' is not installed.", call. = FALSE)
        }
    } else {
        stop("Unknown backend: '", backend, "'.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Run MOFA2 safely with backend switching
#'
#' @param exposomicset_mo A MultiAssayExperiment object
#' @param n_factors Number of latent factors to infer
#' @return A trained MOFA model object
#' @keywords internal
#' @noRd
.run_mofa2 <- function(exposomicset_mo, n_factors) {
    backend <- .mofa_backend()
    message("Using MOFA backend: ", backend)
    .check_mofa_safe(backend)


    # Create MOFA object
    mofa <- MOFA2::create_mofa(exposomicset_mo)

    # Options
    model_opts <- MOFA2::get_default_model_options(mofa)
    model_opts$num_factors <- n_factors
    data_opts <- MOFA2::get_default_data_options(mofa)
    train_opts <- MOFA2::get_default_training_options(mofa)

    # Prepare model
    mofa <- MOFA2::prepare_mofa(
        object = mofa,
        data_options = data_opts,
        model_options = model_opts,
        training_options = train_opts
    )

    outfile <- file.path(tempdir(), "mofa_model.hdf5")

    # Switch runner based on backend
    if (backend == "reticulate") {
        mofa_trained <- MOFA2::run_mofa(mofa, outfile, use_basilisk = FALSE)
    } else {
        mofa_trained <- MOFA2::run_mofa(mofa, outfile, use_basilisk = TRUE)
    }

    # Reload trained model
    MOFA2::load_model(outfile)
}

# --- MCIA Integration Function ----------------
#' @title Run MCIA Integration
#' @description Applies MCIA using `nipalsMCIA` on a
#' MultiAssayExperiment.
#'
#' @param exposomicset_mo A MultiAssayExperiment object containing omics data.
#' @param n_factors Integer. Number of components (PCs) to compute.
#'
#' @return An MCIA result object from `nipalsMCIA::nipals_multiblock()`.
#' @keywords internal
#' @noRd
.run_mcia <- function(
  exposomicset_mo,
  n_factors
) {
    message("Applying MCIA with `nipalsMCIA`")

    # Run NIPALS MCIA on the MultiAssayExperiment
    .check_suggested("nipalsMCIA")
    result <- nipalsMCIA::nipals_multiblock(
        exposomicset_mo,
        col_preproc_method = "colprofile",
        num_PCs = n_factors,
        tol = 1e-12,
        plots = "none"
    )
}
# --- RGCCA Integration Function ---------------
#' @title Run RGCCA Integration
#' @description Applies unsupervised RGCCA to integrate multiple
#' omics blocks.
#'
#' @param exposomicset_mo A MultiAssayExperiment object containing omics data.
#' @param n_factors Integer. Number of components to extract.
#'
#' @return An RGCCA result object from `RGCCA::rgcca()`.
#' @keywords internal
#' @noRd
#' @importFrom RGCCA rgcca
#' @importFrom SummarizedExperiment assay
#' @importFrom MultiAssayExperiment experiments
.run_rgcca <- function(
  exposomicset_mo,
  n_factors
) {
    message("Applying RGCCA integration.")

    # Prepare data: list of assays (samples x features)
    x <- purrr::map(
        MultiAssayExperiment::experiments(exposomicset_mo),
        ~ t(SummarizedExperiment::assay(.x)) # transpose to samples x features
    )

    # Use default shrinkage or manual choice (e.g., 1 for PCA-like)
    tau <- rep(1, length(x)) # no regularization

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
}
# --- DIABLO Integration Function --------------
#' @title Run DIABLO Supervised Integration
#' @description Applies supervised integration using DIABLO (block.splsda)
#' from mixOmics.
#'
#' @param exposomicset_mo A MultiAssayExperiment object containing omics data.
#' @param n_factors Integer. Number of components per block.
#' @param outcome Character. Column name in colData to use as the outcome.
#'
#' @return A DIABLO result object from `mixOmics::block.splsda()`.
#' @keywords internal
#' @noRd
#' @importFrom mixOmics block.splsda
#' @importFrom SummarizedExperiment assay
#' @importFrom MultiAssayExperiment experiments colData
.run_diablo <- function(
  exposomicset_mo,
  n_factors,
  outcome,
  max.iter,
  near.zero.var
) {
    message("Applying DIABLO supervised integration.")

    if (is.null(outcome)) {
        stop("DIABLO requires an outcome variable for supervision.")
    }
    # Ensure that only the common samples are taken
    common <- purrr::map(
        MultiAssayExperiment::experiments(exposomicset_mo),
        ~ colnames(.x)
    ) |>
        (\(lst) Reduce(intersect, lst))()

    # Subset the samples
    exposomicset_mo <- exposomicset_mo[, common]

    # Extract colData and outcome
    y <- MultiAssayExperiment::colData(exposomicset_mo)[[outcome]]
    if (!is.factor(y)) y <- as.factor(y)

    # Prepare blocks (assays) as list of matrices (samples x features)
    blocks <- purrr::map(
        MultiAssayExperiment::experiments(exposomicset_mo),
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
        design = design,
        max.iter = max.iter,
        near.zero.var = near.zero.var
    )
}
