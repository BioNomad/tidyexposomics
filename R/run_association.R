#' Run Association Analysis
#'
#' Perform GLM-based association testing between a specified outcome and features from exposures, omics,
#' latent factors, or GO PCs. Automatically adjusts for covariates and supports both Gaussian and binomial models.
#'
#' @param expomicset A `MultiAssayExperiment` object containing data and metadata.
#' @param outcome The outcome variable name (must be in `colData`).
#' @param source Source of features to test. One of `"omics"`, `"exposures"`, `"factors"`, `"go_pcs"`.
#' @param covariates Optional vector of covariate names to include in the model.
#' @param feature_set Optional character vector of exposure or GO terms to test.
#' @param top_n Optional integer: if using omics source, select top `n` most variable features.
#' @param family GLM family; `"gaussian"` or `"binomial"`.
#' @param correction_method Method for p-value adjustment (default: `"fdr"`).
#' @param action If `"add"` (default), saves results to metadata; else returns results as list.
#' @param min_genes Minimum number of genes required to compute GO PCs.
#' @param feature_col If using GO PCs, the column in `rowData` for matching gene symbols or IDs.
#' @param mirna_assays Optional character vector of assays to exclude when extracting GO terms.
#' @param print Optional variable name to print results for (useful for debugging).
#'
#' @return If `action = "add"`, returns updated `MultiAssayExperiment`. Otherwise, returns a list of:
#' - `results_df`: tidy summary of associations
#' - `covariates`: the covariates used
#' - `model_data`: model matrix used in the GLMs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_association(expomicset, outcome = "BMI", source = "exposures", covariates = c("age", "sex"))
#' }
run_association <- function(
    expomicset,
    outcome,
    source = c("omics", "exposures", "factors", "go_pcs"),
    covariates = NULL,
    feature_set = NULL,
    top_n = NULL,
    family = "gaussian",
    correction_method = "fdr",
    action = "add",
    min_genes = 10,
    feature_col = NULL,
    mirna_assays = NULL,
    print = NULL
) {
  # Validate inputs
  source <- match.arg(source)

  # grab coldata
  data <- expomicset |>
    MultiAssayExperiment::colData() |>
    as.data.frame()

  # scale numeric columns
  data <- data |>
    dplyr::mutate_if(is.numeric, ~ as.numeric(scale(.)))

  # switch based on input
  features_df <- switch(
    source,
    omics = .extract_omics_features(expomicset,
                                    top_n),

    exposures = .extract_exposures(data,
                                   feature_set),

    factors = .extract_latent_factors(expomicset),

    go_pcs = .extract_go_pcs(expomicset,
                             geneset = feature_set,
                             covariates,
                             min_genes = min_genes,
                             feature_col = feature_col,
                             mirna_assays = mirna_assays)
  )

  # create the model data
  if (source == "exposures") {
    model_data <- data
    feature_cols <- feature_set
  } else {
    model_data <- dplyr::left_join(
      tibble::rownames_to_column(data, "id"),
      features_df,
      by = "id"
    ) |>
      tibble::column_to_rownames("id")

    feature_cols <- setdiff(colnames(features_df), "id")
  }


  # ensure outcome variable is factor if binomial
  if (family == "binomial") {

    model_data[[outcome]] <- as.factor(model_data[[outcome]])

    if (length(unique(model_data[[outcome]])) != 2)

      stop("Binary outcome must have exactly 2 levels")
  }

  # ensure that covariates are have more than one value
  if (!is.null(covariates)) {
    covariates <- covariates[covariates %in% colnames(model_data)]
    if (length(covariates) == 0) {
      stop("No valid covariates provided.")
    }
    for (cov in covariates) {
      if (length(unique(model_data[[cov]])) <= 1) {
        stop(paste("Covariate", cov, "must have more than one unique value."))
      }
    }
  }

  message("Running GLMs.")

  feature_cols <- setdiff(colnames(features_df), "id")

  results <- purrr::map_dfr(feature_cols, function(fcol) {
    fmla <- as.formula(
      paste(outcome, "~", fcol,
            if (!is.null(covariates)) paste("+", paste(covariates, collapse = "+")) else "")
    )

    model <- tryCatch(
      glm(fmla, data = model_data, family = family),
      error = function(e) NULL
    )

    # Only continue if model is valid
    if (is.null(model)) return(NULL)

    # extract results
    model_summary <- broom::tidy(model)
    model_filtered <- model_summary |>
      dplyr::filter(term == fcol)

    if (nrow(model_filtered) == 0) {
      message(paste("Skipping", fcol, "- term not found in model (possibly dropped)"))
      return(NULL)
    }

    if(source %in% c("exposures", "factors")){

    }
    # Compute R^2
    if (family == "gaussian") {
      r2 <- summary(model)$r.squared
      adj_r2 <- summary(model)$adj.r.squared

      # Compute R-squared and adjusted R-squared
      n <- nrow(model_data)
      p <- length(coef(model))  # includes intercept
      r2 <- 1 - model$deviance / model$null.deviance
      adj_r2 <- 1 - ((n - 1) / (n - p)) * (1 - r2)

    } else if (family == "binomial") {
      # Compute R-squared and adjusted R-squared
      n <- nrow(model_data)
      p <- length(coef(model))  # includes intercept

      r2 <- 1 - model$deviance / model$null.deviance
      adj_r2 <- 1 - ((n - 1) / (n - p)) * (1 - r2)

    }

    model_filtered$r2 <- r2
    model_filtered$adj_r2 <- adj_r2
    model_filtered
  })


  # adjust p-values
  results <- results |>
    dplyr::mutate(
      p_adjust = p.adjust(
        p.value,
        method = correction_method)) |>
    dplyr::mutate(outcome = outcome)

  if(source %in% c("omics", "factors")){
    # Grab experiment names
    exp_names <- names(MultiAssayExperiment::experiments(expomicset)) |>
      (\(chr) gsub(" ", "_", chr))()

    matched <- data.frame(exp_name=exp_names,
                          exp_name_clean=names(MultiAssayExperiment::experiments(expomicset)))

    results <- results |>
      dplyr::mutate(
        exp_name = stringr::str_extract(term, paste0("^(", paste(exp_names, collapse = "|"), ")")),
        term = dplyr::case_when(
          grepl(paste(exp_names, collapse = "|"), term) ~
            gsub(paste0("(", paste0(exp_names, "_", collapse = "|"), ")"), "", term),
          .default = term
        )
      ) |>
      inner_join(matched,
                 by=c("exp_name"="exp_name")) |>
      dplyr::select(-exp_name) |>
      dplyr::rename(category=exp_name_clean)

    if (source == "factors"){
      if (MultiAssayExperiment::metadata(expomicset) |>
          purrr::pluck("multiomics_integration",
                       "integration_results",
                       "method") %in% c("DIABLO","RGCCA")){
        results <- results |>
          dplyr::mutate(term=paste(category,term,sep=" "))
      }
    }
  }

  if (source == "exposures"){
    results <- results |>
      dplyr::left_join(
        expomicset |>
          MultiAssayExperiment::metadata() |>
          purrr::pluck("codebook") ,
        by=c("term"="variable")
      )
  }

  if (source == "omics"){

    feature_df <- lapply(
      names(
        MultiAssayExperiment::experiments(expom_1)),
      function(name){
        SummarizedExperiment::rowData(
          MultiAssayExperiment::experiments(expom_1)[[name]]) |>
          as.data.frame() |>
          tibble::rownames_to_column(".feature") |>
          mutate(exp_name=name)}) |>
      bind_rows()

    results <- results |>
      dplyr::mutate(
        category = gsub("_", " ", category)
      ) |>
      left_join(
        feature_df,
        by = c("term" = ".feature",
               "category" = "exp_name")
      )
  }

  if (
    (paste0("assoc_", source) %in%
     names(MultiAssayExperiment::metadata(expomicset)$association)) &&
    identical(covariates,MultiAssayExperiment::metadata(expomicset)$association[[paste0("assoc_", source)]]$covariates )){
    if(any(results$term %in%
           MultiAssayExperiment::metadata(expomicset)$association[[paste0("assoc_", source)]]$results_df$term)){
      stop("Association results for this feature already exist in the metadata.")
    } else{
      results <- expomicset |>
        MultiAssayExperiment::metadata() |>
        purrr::pluck("association") |>
        purrr::pluck(paste0("assoc_", source)) |>
        purrr::pluck("results_df") |>
        bind_rows(results)
    }
  }

  if (!is.null(print)) {
    print(results |>
            dplyr::select(term, p.value))
  }

  # add results to metadata
  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$association[[
      paste0("assoc_", source)]] <- list(
        results_df = results,
        covariates = covariates
      )

    # Add in a step record
    step_record <- list(
      run_association = list(
        timestamp = Sys.time(),
        params = list(
          outcome = outcome,
          source = source,
          covariates = covariates,
          feature_set = feature_set,
          top_n = top_n,
          family = family,
          correction_method = correction_method,
          min_genes = min_genes,
          feature_col = feature_col,
          mirna_assays = mirna_assays
        ),
        notes = paste0("Performed association analysis using source: ", source)
      )
    )

    MultiAssayExperiment::metadata(expomicset)$summary$steps <- c(
      MultiAssayExperiment::metadata(expomicset)$summary$steps,
      step_record
    )

    return(expomicset)

  } else {
    return(list(
      results_df = results,
      covariates = covariates,
      model_data = model_data
    ))
  }
}

#' Extract Omics Features
#'
#' Helper for extracting and optionally selecting top omics features from a MultiAssayExperiment
#' for modeling.
#'
#' @keywords internal
#' @noRd
.extract_omics_features <- function(expomicset,
                                    top_n = NULL) {
  log2_assays <- .log2_multiassay(expomicset)

  selected <- if (!is.null(top_n)) {
    .top_var_multiassay(log2_assays, n = top_n)
  }else {
    lapply(MultiAssayExperiment::experiments(log2_assays), rownames)
  }

  scaled <- .scale_multiassay(log2_assays,
                              log2 = FALSE)

  dfs <- purrr::imap(
    MultiAssayExperiment::experiments(scaled),
    function(se, name) {
      df <- SummarizedExperiment::assay(
        se[selected[[name]], , drop = FALSE]) |>
        t() |>
        as.data.frame()

      names(df) <- paste0(name, "_", names(df)) |>
        (\(chr) gsub(" |-", "_", chr))()

      tibble::rownames_to_column(df, "id")
    })

  purrr::reduce(dfs, dplyr::full_join, by = "id")
}


#' Extract Exposure Features
#'
#' Subset exposure features (from colData) for modeling.
#'
#' @keywords internal
#' @noRd
.extract_exposures <- function(data, feature_set) {
  df <- data[, feature_set, drop = FALSE] |> as.data.frame()
  tibble::rownames_to_column(df, "id")
}

#' Extract Latent Factors
#'
#' Retrieve latent factors (MOFA, MCIA, MCCA, DIABLO, or RGCCA) for modeling from metadata.
#'
#' @keywords internal
#' @noRd
.extract_latent_factors <- function(expomicset) {
  result <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("multiomics_integration") |>
    purrr::pluck("integration_results")

  mat <- switch(
    result$method,

    "MOFA" = MOFA2::get_factors(result$result)[[1]],

    "MCIA" = result$result@global_scores,

    "MCCA" = {
      result$result$sample_scores |>
        purrr::map(~ .x |>
                     as.data.frame() |>
                     tibble::rownames_to_column("sample") |>
                     tidyr::pivot_longer(-sample, names_to = "factor", values_to = "weight")) |>
        dplyr::bind_rows(.id = "exp_name") |>
        dplyr::group_by(factor, sample) |>
        dplyr::reframe(weight = mean(weight), .groups = "drop") |>
        tidyr::pivot_wider(names_from = "factor", values_from = "weight") |>
        tibble::column_to_rownames("sample") |>
        (\(df) {colnames(df) = gsub(" ","_",colnames(df));df})()
    },

    "DIABLO" = {
      result$result$variates |>
        (\(lst) lst[names(lst) != "Y"])() |>
        purrr::imap_dfc( ~ {
          comp_names <- paste(.y, colnames(.x), sep = " ")
          df <- as.data.frame(.x)
          colnames(df) <- comp_names
          df
        }) |>
        tibble::rownames_to_column("id") |>
        tibble::column_to_rownames("id") |>
        (\(df) {colnames(df) = gsub(" ","_",colnames(df));df})()
    },

    "RGCCA" = {
      purrr::imap_dfc(result$result$Y, ~ {
        comp_names <- paste(.y, colnames(.x), sep = " ")
        df <- as.data.frame(.x)
        colnames(df) <- comp_names
        df
      }) |>
        tibble::rownames_to_column("id") |>
        tibble::column_to_rownames("id") |>
        (\(df) {colnames(df) = gsub(" ","_",colnames(df));df})()
    },

    stop("Unsupported integration method.")
  )

  mat <- as.data.frame(mat)
  tibble::rownames_to_column(mat, "id")
}



#' Extract GO Principal Components
#'
#' For each GO group in the enrichment results, compute the first PC of the associated genes.
#'
#' @keywords internal
#' @noRd
.extract_go_pcs <- function(
    expomicset, geneset,
    covariates,
    min_genes = 10,
    feature_col = NULL,
    mirna_assays = NULL) {

  enrich_res <- expomicset |>
    MultiAssayExperiment::metadata() |>
    purrr::pluck("enrichment") |>
    purrr::pluck(geneset)

  if (!is.null(mirna_assays)){
    enrich_res <- enrich_res |>
      dplyr::filter(!exp_name %in% mirna_assays)
  }

  pc_dfs <- purrr::pmap_dfr(

    enrich_res |>
      dplyr::distinct(
        exp_name,
        Cluster,
        go_group),

    function(exp_name,
             Cluster,
             go_group) {

      df <- enrich_res |>
        dplyr::filter(exp_name == !!exp_name,
                      Cluster == !!Cluster,
                      go_group == !!go_group)

      genes <- unique(unlist(stringr::str_split(df$geneID, "/")))

      if (length(genes) < min_genes) return(NULL)

      exp <- .update_assay_colData(expomicset, exp_name)

      if (!is.null(feature_col)) {
        genes <- exp |>
          tidybulk::pivot_transcript() |>
          dplyr::filter(!!rlang::sym(feature_col) %in% genes) |>
          dplyr::pull(.feature)

      }

      assay <- SummarizedExperiment::assay(exp)

      assay <- assay[rownames(assay) %in% genes, , drop = FALSE]

      if (nrow(assay) < 2 || all(apply(assay, 1, var) == 0)) return(NULL)

      pc1 <- prcomp(t(log2(assay + 1)), scale. = TRUE)$x[, 1]

      id <- rownames(pc1)

      tibble::tibble(
        id = id,
        !!paste("PC", exp_name, Cluster, go_group, sep = "/") := pc1)
    }
  )

  return(pc_dfs)
}
