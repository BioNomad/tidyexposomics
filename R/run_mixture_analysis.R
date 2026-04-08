#' Run Mixture Analysis
#'
#' Perform mixture analysis to assess joint effects of multiple exposures
#' on an outcome using quantile g-computation (qgcomp), weighted quantile
#' sum (WQS), or Bayesian kernel machine regression (BKMR).
#'
#' @param exposomicset A `MultiAssayExperiment` object.
#' @param outcome The outcome variable name (must be in `colData`).
#' @param exposures Character vector of exposure variable names.
#'   If `NULL`, uses all numeric variables from the codebook.
#' @param covariates Optional character vector of covariate names.
#' @param method One of `"qgcomp"`, `"wqs"`, or `"bkmr"`.
#' @param family `"gaussian"`, `"binomial"`, `"poisson"`, or `"cox"`.
#' @param n_quantiles Number of quantiles. Default is 4.
#' @param n_boot Number of bootstrap samples for qgcomp/WQS. Default is 200.
#'   For qgcomp, set to 0 for no bootstrap (faster, uses Delta method for CI).
#' @param repeat_holdout Number of repeated holdout validations for WQS.
#'   Default is 1.
#' @param validation Proportion of data for validation in WQS. Default is 0.6.
#' @param direction For WQS: `"positive"`, `"negative"`, or `"both"`.
#' @param lambda Penalization parameter for WQS. Default is 0 (no penalty).
#' @param n_iter Number of iterations for BKMR. Default is 10000.
#' @param varsel Logical. Whether to perform variable selection in BKMR.
#'   Default is `TRUE`.
#' @param degree Polynomial degree for qgcomp MSM. Default is 1 (linear).
#' @param id Optional cluster ID variable name for robust standard errors.
#' @param weights Optional sampling weights variable name.
#' @param seed Random seed for reproducibility.
#' @param parallel Logical. Whether to use parallel processing for bootstrap.
#'   Default is `FALSE`.
#' @param keep_fit Logical. Whether to keep the full model fit object.
#'   Default is `FALSE`. Set to `TRUE` if you need access to the underlying
#'   model for diagnostics or plotting.
#' @param action If `"add"` (default), saves results to metadata.
#'   if `"get"`, returns results as list.
#' @param ... Additional arguments passed to the underlying method.
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`.
#'   Otherwise, returns a list containing:
#'
#' For all methods:
#' - `weights`: tibble of exposure weights
#' - `mixture_effect`: tibble with overall mixture effect estimate
#' - `method`, `outcome`, `exposures`, `covariates`: analysis metadata
#'
#' Method-specific outputs:
#' - qgcomp: `partial_effects` tibble
#' - wqs: `model_summary` tibble
#' - bkmr: `pips` tibble (posterior inclusion probabilities)
#'
#' @details
#' Methods available:
#' - `qgcomp`: Quantile g-computation via `qgcomp` package. This is a fast
#'   method and allows negative weights, estimates effects in both directions
#'   simultaneously. Recommended as the default method. Supports clustering via
#'    `id` parameter and sampling weights via `weights` parameter.
#' - `wqs`: Weighted Quantile Sum regression via `gWQS` package. Constrained
#'   weights (need to be all positive or all negative),
#'   this is good for directional hypotheses.
#' - `bkmr`: Bayesian Kernel Machine Regression via `bkmr` package. This is a
#'   flexible nonparametric approach, which can detect nonlinearity and
#'    interactions, but is a slower method.
#'
#' @examples
#' mae <- make_example_data(n_samples = 100, return_mae = TRUE)
#'
#' # Quantile g-computation (fast, recommended)
#' mae |>
#'     run_mixture_analysis(
#'         outcome = "outcome_bmi",
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         method = "qgcomp"
#'     )
#'
#' # qgcomp without bootstrap (even faster)
#' mae |>
#'     run_mixture_analysis(
#'         outcome = "outcome_bmi",
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         method = "qgcomp",
#'         n_boot = 0
#'     )
#'
#' # WQS with both directions
#' mae |>
#'     run_mixture_analysis(
#'         outcome = "outcome_bmi",
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         method = "wqs",
#'         direction = "both"
#'     )
#'
#' @export
run_mixture_analysis <- function(
  exposomicset,
  outcome,
  exposures = NULL,
  covariates = NULL,
  method = c("qgcomp", "wqs", "bkmr"),
  family = "gaussian",
  n_quantiles = 4,
  n_boot = 200,
  repeat_holdout = 1,
  validation = 0.6,
  direction = "both",
  lambda = 0,
  n_iter = 10000,
  varsel = TRUE,
  degree = 1,
  id = NULL,
  weights = NULL,
  seed = 123,
  parallel = FALSE,
  keep_fit = FALSE,
  action = "add",
  ...
) {
    method <- match.arg(method)

    # Check suggested packages
    pkg <- switch(method,
        qgcomp = "qgcomp",
        wqs = "gWQS",
        bkmr = "bkmr"
    )
    .check_suggested(pkg)

    # Get exposures from codebook if not provided
    if (is.null(exposures)) {
        exposures <- MultiAssayExperiment::metadata(exposomicset)$codebook$variable
    }

    # Remove covariates from exposures
    exposures <- setdiff(exposures, covariates)

    # Prepare model data
    model_data <- exposomicset |>
        MultiAssayExperiment::colData() |>
        as.data.frame()

    # Columns needed for complete cases
    vars_needed <- c(outcome, exposures, covariates)
    if (!is.null(id)) vars_needed <- c(vars_needed, id)
    if (!is.null(weights)) vars_needed <- c(vars_needed, weights)

    # Validate
    missing <- setdiff(vars_needed, colnames(model_data))
    if (length(missing) > 0) {
        stop("Variables not found in colData: ", paste(missing, collapse = ", "))
    }

    # Keep only numeric exposures for mixture methods
    numeric_exposures <- exposures[sapply(model_data[exposures], is.numeric)]
    if (length(numeric_exposures) < length(exposures)) {
        dropped <- setdiff(exposures, numeric_exposures)
        message("Dropping non-numeric exposures: ", paste(dropped, collapse = ", "))
        exposures <- numeric_exposures
    }

    # Explicitly create complete case dataset
    # recommended by qgcomp docs
    model_data <- model_data[complete.cases(model_data[, vars_needed]), vars_needed]

    message(sprintf(
        "Running %s with %d exposures, %d samples",
        toupper(method), length(exposures), nrow(model_data)
    ))

    # Run method
    set.seed(seed)
    results <- switch(method,
        qgcomp = .run_qgcomp(
            data = model_data,
            outcome = outcome,
            exposures = exposures,
            covariates = covariates,
            family = family,
            n_quantiles = n_quantiles,
            n_boot = n_boot,
            degree = degree,
            id = id,
            weights = weights,
            seed = seed,
            parallel = parallel,
            keep_fit = keep_fit,
            ...
        ),
        wqs = .run_wqs(
            data = model_data,
            outcome = outcome,
            exposures = exposures,
            covariates = covariates,
            family = family,
            n_quantiles = n_quantiles,
            n_boot = n_boot,
            repeat_holdout = repeat_holdout,
            validation = validation,
            direction = direction,
            lambda = lambda,
            keep_fit = keep_fit,
            ...
        ),
        bkmr = .run_bkmr(
            data = model_data,
            outcome = outcome,
            exposures = exposures,
            covariates = covariates,
            family = family,
            n_iter = n_iter,
            varsel = varsel,
            keep_fit = keep_fit,
            ...
        )
    )

    # Add metadata to results
    results$method <- method
    results$outcome <- outcome
    results$exposures <- exposures
    results$covariates <- covariates
    results$family <- family
    results$n_samples <- nrow(model_data)

    message("Mixture analysis completed.")

    if (action == "add") {
        MultiAssayExperiment::metadata(exposomicset)$mixture_analysis[[method]] <- results

        step_record <- list(
            run_mixture_analysis = list(
                timestamp = Sys.time(),
                params = list(
                    outcome = outcome,
                    exposures = exposures,
                    covariates = covariates,
                    method = method,
                    family = family,
                    n_quantiles = n_quantiles,
                    n_boot = n_boot,
                    seed = seed
                ),
                notes = sprintf("Performed %s mixture analysis", toupper(method))
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )

        return(exposomicset)
    } else {
        return(results)
    }
}

# ---- Quantile G-Computation ----

#' Run quantile g-computation
#' @keywords internal
#' @noRd
.run_qgcomp <- function(
  data,
  outcome,
  exposures,
  covariates,
  family,
  n_quantiles,
  n_boot,
  degree,
  id,
  weights,
  seed,
  parallel,
  keep_fit,
  ...
) {
    .check_suggested("qgcomp")
    suppressPackageStartupMessages(require(qgcomp, quietly = TRUE))

    # Build formula
    if (!is.null(covariates) && length(covariates) > 0) {
        formula <- stats::as.formula(
            paste(outcome, "~", paste(c(exposures, covariates), collapse = " + "))
        )
    } else {
        formula <- stats::as.formula(
            paste(outcome, "~", paste(exposures, collapse = " + "))
        )
    }

    # Extract weights if specified
    wts <- if (!is.null(weights)) data[[weights]] else NULL

    # Choose function based on bootstrap and clustering needs
    if (n_boot > 0) {
        fit <- qgcomp::qgcomp.boot(
            formula,
            data = data,
            expnms = exposures,
            family = family,
            q = n_quantiles,
            B = n_boot,
            degree = degree,
            seed = seed,
            parallel = parallel,
            parplan = parallel,
            ...
        )
    } else if (!is.null(id) || !is.null(wts)) {
        fit <- qgcomp::qgcomp.ee(
            formula,
            data = data,
            expnms = exposures,
            family = family,
            q = n_quantiles,
            id = id,
            weights = wts,
            ...
        )
    } else {
        fit <- qgcomp::qgcomp.noboot(
            formula,
            data = data,
            expnms = exposures,
            family = family,
            q = n_quantiles,
            ...
        )
    }

    .extract_qgcomp_results(fit, exposures, keep_fit)
}

#' Extract qgcomp results
#' @keywords internal
#' @noRd
.extract_qgcomp_results <- function(
  fit,
  exposures,
  keep_fit = FALSE
) {
    # For no boot the weights are stored directly
    # For boot, calculate from coefficients
    if (!is.null(fit$pos.weights)) {
        pos_weights <- fit$pos.weights
        neg_weights <- fit$neg.weights
        pos_psi <- fit$pos.psi
        neg_psi <- fit$neg.psi
        pos_size <- fit$pos.size
        neg_size <- fit$neg.size
    } else {
        # Bootstrap version calculate weights from coefficients
        coefs <- fit$fit$coefficients[exposures]

        pos_coefs <- coefs[coefs > 0]
        pos_sum <- sum(pos_coefs)
        pos_weights <- if (length(pos_coefs) > 0) pos_coefs / pos_sum else numeric(0)
        pos_psi <- pos_sum
        pos_size <- length(pos_coefs)

        neg_coefs <- coefs[coefs < 0]
        neg_sum <- sum(abs(neg_coefs))
        neg_weights <- if (length(neg_coefs) > 0) abs(neg_coefs) / neg_sum else numeric(0)
        neg_psi <- sum(neg_coefs)
        neg_size <- length(neg_coefs)
    }

    # Build weights tibble
    weights <- tibble::tibble(exposure = exposures) |>
        dplyr::mutate(
            weight_pos = dplyr::case_when(
                exposure %in% names(pos_weights) ~ as.numeric(pos_weights[exposure]),
                TRUE ~ 0
            ),
            weight_neg = dplyr::case_when(
                exposure %in% names(neg_weights) ~ as.numeric(neg_weights[exposure]),
                TRUE ~ 0
            ),
            direction = dplyr::case_when(
                weight_pos > 0 ~ "positive",
                weight_neg > 0 ~ "negative",
                TRUE ~ "none"
            )
        )

    # Mixture effect tibble
    if (length(fit$psi) == 1) {
        mixture_effect <- tibble::tibble(
            term = "psi1",
            estimate = fit$psi,
            std_error = sqrt(fit$var.psi),
            lower_ci = fit$ci[1],
            upper_ci = fit$ci[2],
            p_value = fit$pval
        )
    } else {
        mixture_effect <- tibble::tibble(
            term = paste0("psi", seq_along(fit$psi)),
            estimate = fit$psi,
            std_error = sqrt(diag(fit$var.psi)),
            lower_ci = fit$ci[, 1],
            upper_ci = fit$ci[, 2],
            p_value = fit$pval
        )
    }

    # Partial effects tibble
    partial_effects <- tibble::tibble(
        direction = c("positive", "negative"),
        partial_effect = c(pos_psi, neg_psi),
        n_exposures = c(pos_size, neg_size)
    )

    results <- list(
        weights = weights,
        mixture_effect = mixture_effect,
        partial_effects = partial_effects
    )

    if (keep_fit) {
        results$fit <- fit
    }

    results
}

# ---- Weighted Quantile Sum ----

#' Run WQS regression
#' @keywords internal
#' @noRd
.run_wqs <- function(
  data,
  outcome,
  exposures,
  covariates,
  family,
  n_quantiles,
  n_boot,
  repeat_holdout,
  validation,
  direction,
  lambda,
  keep_fit,
  ...
) {
    .check_suggested("gWQS")
    suppressPackageStartupMessages(require(gWQS, quietly = TRUE))

    # Build formula based on direction
    wqs_terms <- switch(direction,
        positive = "wqs",
        negative = "wqs",
        both = "pwqs + nwqs"
    )

    if (!is.null(covariates) && length(covariates) > 0) {
        formula <- stats::as.formula(
            paste(outcome, "~", wqs_terms, "+", paste(covariates, collapse = " + "))
        )
    } else {
        formula <- stats::as.formula(paste(outcome, "~", wqs_terms))
    }

    b1_pos <- direction != "negative"

    fit <- gWQS::gwqs(
        formula = formula,
        mix_name = exposures,
        data = data,
        q = n_quantiles,
        b = n_boot,
        rh = repeat_holdout,
        validation = validation,
        b1_pos = b1_pos,
        family = family,
        lambda = lambda,
        ...
    )

    .extract_wqs_results(fit, direction, keep_fit)
}

#' Extract WQS results
#' @keywords internal
#' @noRd
.extract_wqs_results <- function(
  fit,
  direction,
  keep_fit = FALSE
) {
    # Weights tibble
    weights_df <- fit$final_weights |>
        tibble::as_tibble()

    if (direction == "both") {
        weights <- weights_df |>
            dplyr::transmute(
                exposure = mix_name,
                weight_pos = mean_weight_p,
                weight_neg = mean_weight_n
            )
    } else {
        weights <- weights_df |>
            dplyr::transmute(
                exposure = mix_name,
                weight = mean_weight
            )
    }

    # Model summary tibble
    model_summary <- fit$fit |>
        broom::tidy()

    # Mixture effect tibble for the WQS coefficients
    mixture_effect <- model_summary |>
        dplyr::filter(grepl("wqs", term)) |>
        dplyr::transmute(
            term = term,
            estimate = estimate,
            std_error = std.error,
            statistic = statistic,
            p_value = p.value
        )

    results <- list(
        weights = weights,
        mixture_effect = mixture_effect,
        model_summary = model_summary,
        aic = fit$fit$aic,
        direction = direction
    )

    if (keep_fit) {
        results$fit <- fit
    }

    results
}

# ---- Bayesian Kernel Machine Regression ----

#' Run BKMR
#' @keywords internal
#' @noRd
.run_bkmr <- function(
  data,
  outcome,
  exposures,
  covariates,
  family,
  n_iter,
  varsel,
  keep_fit,
  ...
) {
    .check_suggested("bkmr")

    Z <- as.matrix(data[, exposures, drop = FALSE])
    colnames(Z) <- exposures
    y <- data[[outcome]]

    X <- if (!is.null(covariates) && length(covariates) > 0) {
        as.matrix(data[, covariates, drop = FALSE])
    } else {
        NULL
    }

    fit <- bkmr::kmbayes(
        y = y,
        Z = Z,
        X = X,
        iter = n_iter,
        family = family,
        varsel = varsel,
        verbose = TRUE,
        ...
    )

    results <- list()

    # PIPs tibble (only with variable selection)
    if (varsel) {
        pips_raw <- bkmr::ExtractPIPs(fit)
        results$weights <- tibble::tibble(
            exposure = exposures,
            pip = pips_raw$PIP
        ) |>
            dplyr::arrange(dplyr::desc(pip))
    } else {
        # If no variable selection create weights tibble with NA
        results$weights <- tibble::tibble(
            exposure = exposures,
            pip = NA_real_
        )
    }

    # Posterior summary tibble
    post_summary <- summary(fit)
    if (!is.null(post_summary$summary.betas) && nrow(post_summary$summary.betas) > 0) {
        results$mixture_effect <- post_summary$summary.betas |>
            as.data.frame() |>
            tibble::rownames_to_column("term") |>
            tibble::as_tibble()
    } else {
        results$mixture_effect <- tibble::tibble(
            term = character(),
            estimate = numeric()
        )
    }

    # Include in output if keep_fit = TRUE
    if (keep_fit) {
        results$fit <- fit
    }

    results
}
