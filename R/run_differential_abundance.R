#' Run Differential Abundance Analysis
#'
#' Performs differential abundance testing across all assays in a
#' `MultiAssayExperiment` object using a specified statistical method.
#' The function updates each assay with its
#' corresponding `colData`, fits the model using the provided formula,
#' and combines the results into a unified table.
#'
#' @param exposomicset A `MultiAssayExperiment` containing assays to test.
#' @param formula A model formula for the differential analysis
#' (e.g., ~ group + batch).
#' @param abundance_col Character. The name of the assay matrix to use
#' for abundance values. Default is `"counts"`.
#' @param method Character. Differential analysis method to use.
#' Currently supports `"limma_trend"` (default).
#' @param contrasts A named list of contrasts for pairwise comparisons.
#'  Default is `NULL` (uses default group comparisons).
#' @param scaling_method Character. Scaling method to apply before modeling.
#' Options include `"none"` (default), `"zscore"`, etc.
#' @param action Character. Whether to `"add"` results to `exposomicset`
#'  metadata or `"get"` the results as a data frame. Default is `"add"`.
#'
#' @return Either the updated `MultiAssayExperiment` (if `action = "add"`)
#'  or a tibble with differential abundance results (if `action = "get"`).
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # perform differential abundance analysis
#' mae <- run_differential_abundance(
#'     exposomicset = mae,
#'     formula = ~ smoker + sex,
#'     abundance_col = "counts",
#'     method = "limma_trend",
#'     action = "add"
#' )
#'
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom stats binomial
#' @importFrom rlang .data
#' @export
run_differential_abundance <- function(
    exposomicset,
    formula,
    abundance_col = "counts",
    method = "limma_trend",
    contrasts = NULL,
    scaling_method = "none",
    action = "add") {
    message("Running differential abundance testing.")
    .check_suggested(pkg = "edgeR")
    .check_suggested(pkg = "limma")

    # Initialize a data frame to store results
    da_results_df <- list()

    # Iterate through assays in exposomicset
    for (exp_name in names(MultiAssayExperiment::experiments(exposomicset))) {
        message("Processing assay: ", exp_name)

        # Update assay with colData
        exp <- .update_assay_colData(exposomicset, exp_name)

        # Run differential analysis using `.run_se_differential_abundance`
        res <- .run_se_differential_abundance(
            se = exp,
            formula = formula,
            abundance_col = abundance_col,
            method = method,
            scaling_method = scaling_method,
            contrasts = contrasts
        )

        # If results exist, append assay name and store them
        if (!is.null(res) && nrow(res) > 0) {
            res <- res |>
                dplyr::mutate(exp_name = exp_name)
            da_results_df[[exp_name]] <- res
        } else {
            warning("No significant results found for assay: ", exp_name)
        }
    }

    # Combine results across assays
    final_results <- da_results_df |>
        dplyr::bind_rows() |>
        as_tibble()

    message("Differential abundance testing completed.")

    if (action == "add") {
        all_metadata <- MultiAssayExperiment::metadata(exposomicset)
        all_metadata$differential_analysis$differential_abundance <- final_results
        MultiAssayExperiment::metadata(exposomicset) <- all_metadata

        # Add step record to metadata
        step_record <- list(
            run_differential_abundance = list(
                timestamp = Sys.time(),
                params = list(
                    formula = format(formula),
                    method = method,
                    scaling_method = scaling_method,
                    abundance_col = abundance_col,
                    contrasts = contrasts
                ),
                notes = "Performed differential abundance analysis across all assays."
            )
        )

        MultiAssayExperiment::metadata(exposomicset)$summary$steps <- c(
            MultiAssayExperiment::metadata(exposomicset)$summary$steps,
            step_record
        )
        return(exposomicset)
    } else if (action == "get") {
        return(final_results)
    } else {
        stop("Invalid action specified. Use 'add' or 'get'.")
    }
    return(exposomicset)
}


# --- Limma Trend Option ----------------------
#' Run Differential Abundance Analysis using limma-trend
#'
#' Performs differential abundance analysis on a `SummarizedExperiment`
#' using the limma-trend method. Limma-trend can be used for preprocessed
#' (e.g. log-transformed) abundance data such as logCPM, log-intensity,
#' or log2(TPM+1).
#'
#' @importFrom stats binomial
#' @importFrom rlang .data
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importMethodsFrom tidybulk pivot_transcript
#' @keywords internal
#' @noRd
.run_limma_trend <- function(
    se,
    formula,
    abundance_col,
    contrasts = NULL,
    scaling_method = "none",
    robust = TRUE) {
    .check_suggested("limma")

    mat <- SummarizedExperiment::assay(se, abundance_col)
    if (!is.matrix(mat)) mat <- as.matrix(mat)

    # Optional scaling
    if (identical(scaling_method, "zscore")) {
        mat <- t(scale(t(mat)))
    }

    # Log2 transform if needed
    if (max(mat, na.rm = TRUE) > 50) {
        mat <- log2(mat + 1)
    }

    df <- as.data.frame(SummarizedExperiment::colData(se))
    design <- stats::model.matrix(formula, data = df)

    fit <- limma::lmFit(mat, design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = robust)

    results <- list()
    if (!is.null(contrasts)) {
        cm <- limma::makeContrasts(contrasts = contrasts, levels = design)
        fit2 <- limma::contrasts.fit(fit, cm)
        fit2 <- limma::eBayes(fit2, trend = TRUE, robust = robust)
        for (i in seq_len(ncol(cm))) {
            tb <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none")
            tb$contrast <- colnames(cm)[i]
            results[[i]] <- tb
        }
    } else {
        coef_idx <- setdiff(
            seq_len(ncol(design)),
            which(colnames(design) == "(Intercept)")
        )
        for (k in coef_idx) {
            tb <- limma::topTable(fit, coef = k, number = Inf, sort.by = "none")
            tb$contrast <- colnames(design)[k]
            results[[length(results) + 1]] <- tb
        }
    }

    out <- dplyr::bind_rows(results)

    out <- out |>
        tibble::rownames_to_column("feature") |>
        dplyr::mutate(
            feature = gsub("\\.\\.\\..*", "", feature),
            method = "limma_trend",
            .abundant = TRUE,
            scaling = scaling_method
        ) |>
        dplyr::left_join(
            se |>
                pivot_transcript(),
            by = c("feature" = ".feature")
        ) |>
        dplyr::filter(grepl(all.vars(formula)[1], contrast)) |>
        tibble::as_tibble()

    return(out)
}
