# --- Define Custom Color Pallette  ---

tidy_exp_pal <- c(
  "#1E0C1F",
  "#00c9c1",
  "#cf68d9",
  "#6d7ea9",
  "#e8d5b5",
  "#944470",
  "#67abff",
  "#b5a7b6",
  "#002f4f",
  "#f2b459",
  "#4f4350",
  "#51aabb",
  "#994a4f",
  "#00a3cd",
  "#240d00",
  "#48817d"
)

divergent_cols <- c('#0000FF','#FF0000','#00FF00','#000033','#FF00B6','#005300','#FFD300','#009FFF','#9A4D42','#00FFBE','#783FC1','#1F9698','#FFACFD','#B1CC71','#F1085C','#FE8F42','#DD00FF','#201A01','#720055','#766C95','#02AD24','#C8FF00','#886C00','#FFB79F','#858567','#A10300','#14F9FF','#00479E','#DC5E93','#93D4FF','#004CFF')

scale_fill_tidy_exp <- function(...,rev = F) {
  if (rev) {
    tidy_exp_pal <- rev(tidy_exp_pal)
  }
  ggplot2::scale_fill_manual(values = tidy_exp_pal, ...)
}

scale_color_tidy_exp <- function(..., rev = F) {
  if (rev) {
    tidy_exp_pal <- rev(tidy_exp_pal)
  }
  ggplot2::scale_color_manual(values = tidy_exp_pal, ...)
}


# --- Update Assay ColData ------

#' Update Assay colData in a MultiAssayExperiment
#'
#' Synchronizes the sample metadata (`colData`) of a specific assay with the global `colData`
#' from a `MultiAssayExperiment` object.
#'
#' @keywords internal
#' @noRd
.update_assay_colData <- function(expomicset,
                                  exp_name) {

  # Retrieve the assay
  assay <- MultiAssayExperiment::experiments(expomicset)[[exp_name]]

  # Extract colData for the assay's samples
  assay_samples <- colnames(assay)
  global_coldata <- as.data.frame(MultiAssayExperiment::colData(expomicset))
  coldata <- global_coldata[rownames(global_coldata) %in% assay_samples, , drop = FALSE]

  # Ensure the sample order matches
  coldata <- coldata[match(assay_samples, rownames(coldata)), , drop = FALSE]

  # Add a check to ensure the order is correct
  if (!identical(rownames(coldata), assay_samples)) {
    stop("Sample order mismatch detected in assay: ", exp_name,
         "\nEnsure the samples in colData are aligned with the assay samples.")
  }

  # Update colData in the assay
  MultiAssayExperiment::colData(assay) <- S4Vectors::DataFrame(coldata)

  return(assay)
}
# --- Log2 MultiAssayExperiment Assays --------
#' Log2 Transform Assays in a MultiAssayExperiment
#'
#' Applies log2 transformation to each assay in a `MultiAssayExperiment` object.
#'
#' @keywords internal
#' @noRd
.log2_multiassay <- function(expomicset) {

  message("Log2-Transforming each assay in MultiAssayExperiment.")

  # Apply log2 transformation to each assay
  log2_experiments <- lapply(MultiAssayExperiment::experiments(expomicset),
                             function(assay_obj) {
                               if (inherits(assay_obj, "SummarizedExperiment")) {
                                 assay_mat <- SummarizedExperiment::assay(assay_obj)
                                 log2_mat <- log2(assay_mat + 1)  # Log2 transformation
                                 SummarizedExperiment::assay(assay_obj) <- log2_mat
                                 return(assay_obj)
                               } else if (is.matrix(assay_obj)) {
                                 return(log2(assay_obj + 1))  # Directly log2 transform matrices
                               } else {
                                 stop("Unsupported assay type. Only SummarizedExperiment and matrices are supported.")
                               }
                             })

  log2_expomicset <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = log2_experiments,
    colData = MultiAssayExperiment::colData(expomicset),
    metadata = MultiAssayExperiment::metadata(expomicset))


  return(log2_expomicset)
}
# --- Scale MultiAssayExperiment Assays --------
#' Scale Assays in a MultiAssayExperiment
#'
#' Standardizes each assay in a `MultiAssayExperiment` object using Z-score scaling.
#'
#' @keywords internal
#' @noRd
.scale_multiassay <- function(expomicset,
                              log2=FALSE) {

  message("Scaling each assay in MultiAssayExperiment.")

  # Check if log2 transformation is needed
  if (log2){
    expomicset <- .log2_multiassay(expomicset)
  }

  # Apply scaling to each assay
  scaled_experiments <- lapply(MultiAssayExperiment::experiments(expomicset),
                               function(assay_obj) {
    if (inherits(assay_obj, "SummarizedExperiment")) {
      assay_mat <- SummarizedExperiment::assay(assay_obj)
      scaled_mat <- scale(assay_mat)  # Standardize (Z-score)
      SummarizedExperiment::assay(assay_obj) <- scaled_mat
      return(assay_obj)
    } else if (is.matrix(assay_obj)) {
      return(scale(assay_obj))  # Directly scale matrices
    } else {
      stop("Unsupported assay type. Only SummarizedExperiment and matrices are supported.")
    }
  })

  # Create a new MultiAssayExperiment with scaled data
  scaled_expomicset <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = scaled_experiments,
    colData = MultiAssayExperiment::colData(expomicset),
    metadata = MultiAssayExperiment::metadata(expomicset))

  return(scaled_expomicset)
}
# --- Get Top Variant Features ---------------
#' Get Top Variant Features
#'
#' Extracts the top features based on variance from a `MultiAssayExperiment` object.
#'
#' @keywords internal
#' @noRd
.top_var_multiassay <- function(expomicset,
                                  n = 1000,
                                  assay_name = NULL) {

  #expomicset <- .log2_multiassay(expomicset)

  # Grab the top n features based on variance
  top_var_experiments <- lapply(MultiAssayExperiment::experiments(expomicset),
                               function(assay_obj) {
                                 if (inherits(assay_obj, "SummarizedExperiment")) {
                                   assay_mat <- SummarizedExperiment::assay(assay_obj)
                                   feature_variance <- apply(assay_mat, 1, var) |>
                                     sort() |>
                                     tail(n=n) |>
                                     names()


                                   return(feature_variance)
                                 } else if (is.matrix(assay_obj)) {

                                   return(apply(assay_obj, 1, var) |>
                                            sort() |>
                                            tail(n=n) |>
                                            names())  # Directly scale matrices
                                 } else {
                                   stop("Unsupported assay type. Only SummarizedExperiment and matrices are supported.")
                                 }
                               })


  return(top_var_experiments)
}
# --- Run Differential Abundance Analysis ------
#' Run Differential Abundance Analysis on a SummarizedExperiment
#'
#' Performs differential abundance analysis using `tidybulk`, with optional contrast testing.
#'
#' @keywords internal
#' @noRd
.run_se_differential_abundance <- function(
    se,
    formula,
    abundance_col = "counts",
    method = "limma_voom",
    scaling_method = "none",
    contrasts = NULL
) {

  # Confirm there are no negative values
  if (min(SummarizedExperiment::assay(se,abundance_col)) < 0) {
    # add a absolute value of minimum and add pseudocount to avoid negative values
    message("Negative values detected. Adding pseudocount to all values..")
    min_value <- abs(min(SummarizedExperiment::assay(se,abundance_col))) + 1
    SummarizedExperiment::assay(se,abundance_col) <- SummarizedExperiment::assay(se,abundance_col) + min_value

  }

  # Check for contrast input
  if (!is.null(contrasts)) {
    res_list <- list()
    for (contrast in contrasts) {

      # Run differential abundance analysis
      contrast_results <- se |>
        # Hard coding minimum counts and proportion to 0 since filtering is per omic
        tidybulk::identify_abundant(minimum_counts = 0,
                          minimum_proportion = 0) |>
        tidybulk::test_differential_abundance(
          formula,
          .abundance = !!sym(abundance_col),
          method = method,
          contrasts = contrast,
          scaling_method = scaling_method
        )

      # Extract results
      res <- as.data.frame(S4Vectors::elementMetadata(contrast_results))
      colnames(res) <- gsub("__.*", "", colnames(res))

      # Add metadata
      res <- res |>
        dplyr::mutate(
          feature = rownames(contrast_results),
          contrast = contrast,
          method = method,
          scaling = scaling_method
        )

      res_list[[contrast]] <- res
    }
    return(dplyr::bind_rows(res_list))
  } else {
    contrast_results <- se |>
      # Hard coding minimum counts and proportion to 0 since filtering is per omic
      tidybulk::identify_abundant(
        minimum_counts = 0,
        minimum_proportion = 0) |>
      tidybulk::test_differential_abundance(
        formula,
        .abundance = !!sym(abundance_col),
        method = method,
        scaling_method = scaling_method
      )

    # Extract results
    res <- as.data.frame(S4Vectors::elementMetadata(contrast_results))
    colnames(res) <- gsub("__.*", "", colnames(res))

    # Add metadata
    res <- res |>
      dplyr::mutate(
        feature = rownames(contrast_results),
        contrast = all.vars(formula)[1],
        method = method,
        scaling = scaling_method
      )

    return(res)
  }
}

# --- Calculate Feature Stability -------
#' Calculate Feature Stability Across Sensitivity Conditions
#'
#' Computes a stability score for features based on their significance across multiple sensitivity tests.
#'
#' @keywords internal
#' @noRd
.calculate_feature_stability <- function(sensitivity_df,
                                         pval_col = "adj.P.Val",
                                         logfc_col = "logFC",
                                         pval_threshold = 0.05) {
  # Extract sensitivity analysis results
  sensitivity_df <- sensitivity_df

  message("Computing feature stability across sensitivity conditions.")

  feature_stability_df <- sensitivity_df |>
    group_by(feature, exp_name) |>
    summarise(
      # selection frequency (fraction of significant p-values)
      presence_rate = mean(!!sym(pval_col) < pval_threshold, na.rm = TRUE),

      # effect size consistency (inverse of coefficient of variation)
      effect_consistency = 1 / (1 + (sd(!!sym(logfc_col), na.rm = TRUE) /
                                       mean(abs(!!sym(logfc_col)), na.rm = TRUE))),

      # main hybrid score combining presence and effect stability
      stability_score = presence_rate * effect_consistency,

      # statistical signal strength based on p-values
      mean_log_p = mean(-log10(!!sym(pval_col) + 1e-10), na.rm = TRUE),

      # hybrid score using log p-values
      logp_weighted_score = mean_log_p * effect_consistency,

      # effect size variability
      sd_logFC = sd(!!sym(logfc_col), na.rm = TRUE),
      iqr_logFC = IQR(!!sym(logfc_col), na.rm = TRUE),
      cv_logFC = sd(!!sym(logfc_col), na.rm = TRUE) / abs(mean(!!sym(logfc_col), na.rm = TRUE)),

      # sign flip frequency of logfc across runs
      sign_flip_freq = mean(sign(!!sym(logfc_col)) != sign(mean(!!sym(logfc_col), na.rm = TRUE)), na.rm = TRUE),

      # variability of p-values on the log scale
      sd_log_p = sd(log10(!!sym(pval_col) + 1e-10), na.rm = TRUE)
    )

  message("Feature stability analysis completed.")
  return(feature_stability_df)
}

# --- Correlate Se with colData --------
#' Correlate SummarizedExperiment Features with Exposure Variables
#'
#' Computes correlations between assay features and exposure variables using Spearman or other correlation methods.
#'
#' @keywords internal
#' @noRd
.correlate_se_with_coldata <- function(
    se,
    exposure_cols,
    correlation_method = "spearman",
    correlation_cutoff = 0.0,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05
) {
  require(tidyverse)
  require(Hmisc)
  require(reshape2)

  message("Performing correlation analysis on summarized experiment.")

  # Ensure colData has the specified exposures
  exposure_data <- SummarizedExperiment::colData(se) |>
    as.data.frame()

  numeric_exposures <- intersect(colnames(exposure_data), exposure_cols)

  if (length(numeric_exposures) == 0) {
    stop("No valid exposure variables found in colData.")
  }

  # Extract assay data
  assay_data <- SummarizedExperiment::assays(se)[[1]] |>
    t() |>
    as.data.frame()

  # Ensure colData and assay samples match
  common_samples <- intersect(rownames(assay_data), rownames(exposure_data))

  if (length(common_samples) == 0) {
    stop("No common samples between assay and colData.")
  }

  # Subset both data frames to common samples
  assay_data <- assay_data[common_samples, , drop = FALSE]
  exposure_data <- exposure_data[common_samples, numeric_exposures, drop = FALSE]

  # Merge into a single dataframe for correlation
  merged_data <- assay_data |>
    tibble::rownames_to_column("id") |>
    dplyr::inner_join(
      exposure_data |>
        tibble::rownames_to_column("id"),
      by = "id") |>
    tibble::column_to_rownames("id")

  # Perform correlation analysis
  message("Running Spearman correlation analysis.")
  correlation_matrix <- merged_data |>
    as.matrix() |>
    Hmisc::rcorr(type = correlation_method)

  # Convert correlation and p-values to tidy format
  correlation_df <- correlation_matrix$r |>
    as.data.frame() |>
    tibble::rownames_to_column("feature") |>
    reshape2::melt(id.vars = "feature") |>
    `colnames<-`(c("feature", "exposure", "correlation"))

  pvalue_df <- correlation_matrix$P |>
    as.data.frame() |>
    tibble::rownames_to_column("feature") |>
    reshape2::melt(id.vars = "feature") |>
    `colnames<-`(c("feature", "exposure", "p.value"))

  # Merge correlation results
  correlation_results <- correlation_df |>
    dplyr::inner_join(pvalue_df,
                      by = c("feature", "exposure")) |>
    dplyr::filter(!(abs(correlation) > 1)) |>
    dplyr::filter(abs(correlation) > correlation_cutoff) |>
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
    dplyr::filter(!!sym(cor_pval_column) < pval_cutoff) |>
    dplyr::arrange(desc(abs(correlation))) |>
    dplyr::filter(feature %in% rownames(se)) |>
    dplyr::filter(exposure %in% numeric_exposures)

  message("Correlation analysis completed.")

  return(correlation_results)
}

# --- Get miRNA Targets ----------
#' Retrieve miRNA Target Genes from Omnipath
#'
#' Queries the Omnipath database to extract target genes for specified miRNAs.
#'
#' @keywords internal
#' @noRd
.get_mirna_targets <- function(mirnas) {
  # pull omnipath database
  mirna_db <- OmnipathR::import_mirnatarget_interactions()
  mirna_targets <- mirna_db |>
    dplyr::filter(source_genesymbol %in% mirnas) |>
    dplyr::select(source_genesymbol,target_genesymbol)
  return(mirna_targets)
}

# --- Differential Abundance Exposure Functional Enrichment ------------
#' Perform Functional Enrichment Analysis on Differentially Abundant Features
#'
#' Conducts enrichment analysis for differentially abundant features correlated with exposures,
#' integrating miRNA target mapping and pathway annotation.
#'
#' @keywords internal
#' @noRd
.da_exposure_functional_enrichment <- function(
    expomicset,
    feature_col="feature",
    mirna_assays = NULL,
    pval_col = "adj.P.Val",
    pval_threshold = 0.05,
    logfc_col = "logFC",
    logfc_threshold = log2(1),
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {

  # Check for differential abundance results
  if (!"differential_abundance" %in% names(MultiAssayExperiment::metadata(expomicset)$differential_analysis)) {
    stop("Please run `run_differential_abundance() first.`")
  }

  # Check for correlation results
  if (!"degs" %in% names(MultiAssayExperiment::metadata(expomicset)$correlation)) {
    stop("Please run `run_correlation() first.`")
  }

  # Extract differential abundance and correlation results
  da_res <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance
  da_cor_res <- MultiAssayExperiment::metadata(expomicset)$correlation$degs

  # Filter and merge results based on p-value and logFC thresholds
  da_res_cor_merged <- da_res |>
    dplyr::filter(!!sym(pval_col) < pval_threshold) |>
    dplyr::filter(abs(!!sym(logfc_col)) > logfc_threshold) |>
    dplyr::mutate(direction=ifelse(logFC>0,"up","down")) |>
    dplyr::select(feature,
                  exp_name,
                  !!sym(feature_col),
                  direction) |>
    dplyr::inner_join(da_cor_res,
               by=c("feature",
                    "exp_name")) |>
    dplyr::mutate(feature_col=!!sym(feature_col))

  if (!is.null(mirna_assays)){
    # Retrieve miRNA target genes and merge with differential abundance results
    mirna_df <- da_res_cor_merged |>
      dplyr::filter(exp_name %in% mirna_assays) |>
      (\(df) split(df,df$exp_name))() |>
      purrr::map(~ .get_mirna_targets(
        .x |>
          dplyr::pull(feature_col)) |>
          dplyr::mutate(direction="down")) |>
      dplyr::bind_rows(.id="exp_name") |>
      dplyr::inner_join(da_res_cor_merged |>
                          dplyr::select(-direction),
                        by=c("source_genesymbol"="feature_col",
                             "exp_name")) |>
      dplyr::select(-source_genesymbol) |>
      dplyr::rename(feature_col=target_genesymbol)

    # Combine results with miRNA target df
    da_res_cor_merged_mapped <- da_res_cor_merged |>
      dplyr::filter(!exp_name %in%
                      c(mirna_assays)) |>
      dplyr::bind_rows(mirna_df)
  } else{
    da_res_cor_merged_mapped <- da_res_cor_merged
  }

  # Establish the universe of features for each assay
  universe_per_assay <- as.list(unique(da_res_cor_merged$exp_name)) |>
    purrr::map( ~
           {
           if (.x %in% mirna_assays) {
             # Retrieve features for each assay
             df <- data.frame(all_features=
                                MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                rownames(),
                              exp_name=.x)

             # Get miRNA target genes
             all_targets <- .get_mirna_targets(df$all_features) |>
               dplyr::pull(target_genesymbol)

             # Create a dataframe with target genes and other features
             df <- data.frame(all_features=all_targets,
                              exp_name=.x)
           } else{
             # Extract features for universe based on feature_col
             df <- data.frame(all_features=
                                MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                tidybulk::pivot_transcript() |>
                                dplyr::mutate(feature=.feature) |>
                                dplyr::pull(!!sym(feature_col)),
                              exp_name=.x)
           }
           return(df)}
    ) |>
    dplyr::bind_rows()


  # Perform functional enrichment analysis using compareCluster
  enrich_res <- da_res_cor_merged_mapped |>
    (\(df) split(df,df$exp_name))() |>
    purrr::map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- clusterProfiler::compareCluster(
        feature_col~direction+category,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          dplyr::filter(exp_name==unique(.x$exp_name)) |>
          dplyr::pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })


  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")

  return(enrich_df)
}

# --- Differential Abundance Functional Enrichment ------------
#' Perform Functional Enrichment Analysis on Differentially Abundant Features
#'
#' Conducts pathway enrichment analysis for differentially abundant features,
#' integrating miRNA target mapping when applicable.
#'
#' @keywords internal
#' @noRd
.da_functional_enrichment <- function(
    expomicset,
    robust = T,
    score_col = "stability_score",
    score_threshold = NULL,
    feature_col="feature",
    mirna_assays = NULL,
    pval_col = "adj.P.Val",
    pval_threshold = 0.05,
    logfc_col = "logFC",
    logfc_threshold = log2(1),
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {
  require(clusterProfiler)
  require(biomaRt)
  require(tidyverse)
  require(org.Hs.eg.db)

  # Check for differential abundance results
  if (!"differential_abundance" %in% names(expomicset@metadata$differential_analysis)) {
    stop("Please run `run_differential_abundance() first.`")
  }

  # Filter and merge results based on p-value and logFC thresholds
  da_res <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$differential_abundance |>
    dplyr::filter(!!sym(pval_col) < pval_threshold) |>
    dplyr::filter(abs(!!sym(logfc_col)) > logfc_threshold) |>
    dplyr::mutate(direction=ifelse(logFC>0,"up","down")) |>
    dplyr::select(feature,
                  !!sym(feature_col),
                  exp_name,
                  direction) |>
    dplyr::mutate(feature_col=!!sym(feature_col))

  # If robust analysis is requested, filter based on stability score
  if (robust) {
    if (!"sensitivity_analysis" %in% names(MultiAssayExperiment::metadata(expomicset)$differential_analysis)) {
      stop("Please run `run_sensitivity_analysis() first.`")
    }

    if(is.null(score_threshold)) {
      score_threshold <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis$score_thresh
    }

    # Filter stability scores based on the specified threshold
    stability_df <- MultiAssayExperiment::metadata(expomicset)$differential_analysis$sensitivity_analysis$feature_stability |>
      dplyr::filter(!!dplyr::sym(score_col) > score_threshold) |>
      dplyr::select(feature,
                    exp_name,
                    !!sym(score_col))

    da_res <- da_res |>
      dplyr::inner_join(stability_df,
                        by = c("feature",
                               "exp_name"))
  }

  if (!is.null(mirna_assays)) {
    # Retrieve miRNA target genes and merge with differential abundance results
    mirna_df <- da_res |>
      filter(exp_name %in% mirna_assays) |>
      (\(df) split(df,df$exp_name))() |>
      purrr::map(~ .get_mirna_targets(
        .x |>
          dplyr::pull(feature_col)) |>
          dplyr::mutate(direction="down")) |>
      dplyr::bind_rows(.id="exp_name") |>
      dplyr::inner_join(da_res |>
                          dplyr::select(-direction),
                        by=c("source_genesymbol"="feature",
                             "exp_name")) |>
      dplyr::select(-source_genesymbol) |>
      dplyr::rename(feature_col=target_genesymbol)

    # Combine results with miRNA target df
    da_res_mapped <- da_res |>
      dplyr::filter(!exp_name %in% c(mirna_assays)) |>
      dplyr::bind_rows(mirna_df)
  } else{
    da_res_mapped <- da_res
  }


  universe_per_assay <- as.list(unique(da_res_cor_merged$exp_name)) |>
    purrr::map( ~
                  {
                    if (.x %in% mirna_assays) {
                      # Retrieve features for each assay
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         rownames(),
                                       exp_name=.x)

                      # Get miRNA target genes
                      all_targets <- .get_mirna_targets(df$all_features) |>
                        dplyr::pull(target_genesymbol)

                      # Create a dataframe with target genes and other features
                      df <- data.frame(all_features=all_targets,
                                       exp_name=.x)
                    } else{
                      # Extract features for universe based on feature_col
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         pivot_transcript() |>
                                         mutate(feature=.feature) |>
                                         dplyr::pull(!!sym(feature_col)),
                                       exp_name=.x)
                    }
                    return(df)}
    ) |>
    dplyr::bind_rows()

  # Perform functional enrichment analysis using compareCluster
  enrich_res <- da_res_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- clusterProfiler::compareCluster(
        feature_col~direction,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          dplyr::filter(exp_name==unique(.x$exp_name)) |>
          dplyr::pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })


  enrich_df <- enrich_res |>
    purrr::map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    dplyr::bind_rows(.id="exp_name")

  return(enrich_df)
}

# --- Factor Feature Exposure Functional Enrichment ----------
#' Perform Functional Enrichment Analysis on Factor-Associated Features
#'
#' Conducts enrichment analysis for top factor-contributing features correlated with exposures,
#' integrating miRNA target mapping when applicable.
#'
#' @keywords internal
#' @noRd
.factor_exposure_functional_enrichment <- function(
    expomicset,
    feature_col="feature",
    mirna_assays = NULL,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {

  # Check for top factor features and correlation results
  if (!"top_factor_features" %in% names(MultiAssayExperiment::metadata(expomicset)$multiomics_integration)) {
    stop("Please run `extract_top_factor_features() first.`")
  }

  if (!"factors" %in% names(MultiAssayExperiment::metadata(expomicset)$correlation)) {
    stop("Please run `run_correlation() first.`")
  }

  # Extract factor results and correlation results
  factor_res <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$top_factor_features
  factor_cor_res <- MultiAssayExperiment::metadata(expomicset)$correlation$factors

  # Filter and merge results based on p-value and logFC thresholds
  factor_res_cor_merged <- factor_res |>
    dplyr::select(feature,
                  exp_name) |>
    dplyr::inner_join(factor_cor_res,
               by=c("feature",
                    "exp_name")) |>
    dplyr::inner_join(pivot_feature(expomicset),
               by=c("feature"=".feature",
                    "exp_name"=".exp_name")) |>
    dplyr::mutate(feature_col=!!sym(feature_col))

  if(!is.null(mirna_assays)){
    # Retrieve miRNA target genes and merge with factor results
    mirna_df <- factor_res_cor_merged |>
      dplyr::filter(exp_name %in% mirna_assays) |>
      (\(df) split(df,df$exp_name))() |>
      map(~ .get_mirna_targets(
        .x |>
          dplyr::pull(feature))) |>
      dplyr::bind_rows(.id="exp_name") |>
      dplyr::inner_join(factor_res_cor_merged,
                 by=c("source_genesymbol"="feature",
                      "exp_name")) |>
      dplyr::select(-source_genesymbol) |>
      dplyr::rename(feature=target_genesymbol)

    # Combine results with miRNA target df
    factor_res_cor_merged_mapped <- factor_res_cor_merged |>
      dplyr::filter(!exp_name %in% c(mirna_assays)) |>
      dplyr::bind_rows(mirna_df)
  }else{
    factor_res_cor_merged_mapped <- factor_res_cor_merged
  }

  # Establish the universe of features for each assay
  universe_per_assay <- as.list(unique(factor_res_cor_merged$exp_name)) |>
    purrr::map( ~
                  {
                    if (.x %in% mirna_assays) {
                      # Retrieve features for each assay
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         rownames(),
                                       exp_name=.x)

                      # Get miRNA target genes
                      all_targets <- .get_mirna_targets(df$all_features) |>
                        dplyr::pull(target_genesymbol)

                      # Create a dataframe with target genes and other features
                      df <- data.frame(all_features=all_targets,
                                       exp_name=.x)
                    } else{
                      # Extract features for universe based on feature_col
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         tidybulk::pivot_transcript() |>
                                         dplyr::mutate(feature=.feature) |>
                                         dplyr::pull(!!sym(feature_col)),
                                       exp_name=.x)
                    }
                    return(df)}
    ) |>
    dplyr::bind_rows()

  # Perform functional enrichment analysis using compareCluster
  enrich_res <- factor_res_cor_merged_mapped |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- compareCluster(
        feature_col~category,
        data = .x,
        fun = fun,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          dplyr::filter(exp_name==unique(.x$exp_name)) |>
          dplyr::pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    })


  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@compareClusterResult
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")

  return(enrich_df)
}

# --- Factor Feature Functional Enrichment ------------
#' Perform Functional Enrichment Analysis on Factor-Associated Features
#'
#' Conducts Gene Ontology (GO) enrichment analysis for top factor-contributing features,
#' integrating miRNA target mapping when applicable.
#'
#' @keywords internal
#' @noRd
.factor_functional_enrichment <- function(
    expomicset,
    feature_col="feature",
    mirna_assays = NULL,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP") {

  # Check for top factor features and correlation results
  if (!"top_factor_features" %in% names(MultiAssayExperiment::metadata(expomicset)$multiomics_integration)) {
    stop("Please run `extract_top_factor_features() first.`")
  }

  # Extract factor results
  factor_res <- MultiAssayExperiment::metadata(expomicset)$multiomics_integration$top_factor_features |>
    dplyr::select(feature,
                  exp_name) |>
    dplyr::inner_join(pivot_feature(expomicset),
               by=c("feature"=".feature",
                    "exp_name"=".exp_name")) |>
    dplyr::mutate(feature_col=!!sym(feature_col))

  if(!is.null(mirna_assays)){
    # Retrieve miRNA target genes and merge with factor results
    mirna_df <- factor_res |>
      filter(exp_name %in% mirna_assays) |>
      (\(df) split(df,df$exp_name))() |>
      map(~ .get_mirna_targets(
        .x |>
          dplyr::pull(feature))) |>
      dplyr::bind_rows(.id="exp_name") |>
      dplyr::inner_join(factor_res,
                 by=c("source_genesymbol"="feature",
                      "exp_name")) |>
      dplyr::select(-source_genesymbol) |>
      dplyr::rename(feature=target_genesymbol)

    # Combine results with miRNA target df
    factor_res <- factor_res |>
      dplyr::filter(!exp_name %in% c(mirna_assays)) |>
      dplyr::bind_rows(mirna_df)
  }else{
    factor_res <- factor_res
  }

  # Establish the universe of features for each assay
  universe_per_assay <- as.list(unique(factor_res$exp_name)) |>
    purrr::map( ~
                  {
                    if (.x %in% mirna_assays) {
                      # Retrieve features for each assay
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         rownames(),
                                       exp_name=.x)

                      # Get miRNA target genes
                      all_targets <- .get_mirna_targets(df$all_features) |>
                        dplyr::pull(target_genesymbol)

                      # Create a dataframe with target genes and other features
                      df <- data.frame(all_features=all_targets,
                                       exp_name=.x)
                    } else{
                      df <- data.frame(all_features=
                                         MultiAssayExperiment::experiments(expomicset)[[.x]] |>
                                         tidybulk::pivot_transcript() |>
                                         dplyr::mutate(feature=.feature) |>
                                         dplyr::pull(!!sym(feature_col)),
                                       exp_name=.x)
                    }
                    return(df)}
    ) |>
    dplyr::bind_rows()

  # Perform functional enrichment analysis using compareCluster
  enrich_res <- factor_res |>
    (\(df) split(df,df$exp_name))() |>
    map(~ {
      message("Working on: ",unique(.x$exp_name))
      enrich <- enrichGO(
        gene = .x$feature_col,
        OrgDb = OrgDb,
        keyType = keyType,
        ont = ont,
        universe = universe_per_assay |>
          dplyr::filter(exp_name==unique(.x$exp_name)) |>
          dplyr::pull(all_features),
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff

      )
    })


  enrich_df <- enrich_res |>
    map(~ {
      if(!is.null(.x)){
        enrich <- .x@result
      } else{
        enrich <- NULL
      }
    }) |>
    bind_rows(.id="exp_name")

  return(enrich_df)
}

# --- Pairwise Overlaps -------------
#' Compute Pairwise Overlaps Between Sets
#'
#' Calculates pairwise overlaps, Jaccard indices, and shared elements for a list of unique sets.
#'
#' @keywords internal
#' @noRd
.get_pairwise_overlaps <- function(sets) {
  # credit for most of the code:
  # https://blog.jdblischak.com/posts/pairwise-overlaps/
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x) length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }

  n_sets <- length(sets)
  set_names <- names(sets)
  n_overlaps <- choose(n = n_sets, k = 2)

  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_num_shared <- integer(length = n_overlaps)
  vec_overlap <- numeric(length = n_overlaps)
  vec_jaccard <- numeric(length = n_overlaps)
  vec_shared_terms <- character(length = n_overlaps)
  overlaps_index <- 1

  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]

      shared_terms <- paste(Reduce(intersect,list(set1,set2)),collapse = ",")

      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                    fromLast = FALSE, nmax = NA))
      num_shared <- length(set_intersect)
      overlap <- num_shared / min(length(set1), length(set2))
      jaccard <- num_shared / length(set_union)

      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_num_shared[overlaps_index] <- num_shared
      vec_overlap[overlaps_index] <- overlap
      vec_jaccard[overlaps_index] <- jaccard
      vec_shared_terms[overlaps_index] <- shared_terms

      overlaps_index <- overlaps_index + 1
    }
  }

  result <- data.frame(source = vec_name1,
                       target = vec_name2,
                       num_shared = vec_num_shared,
                       overlap = vec_overlap,
                       jaccard = vec_jaccard,
                       shared_terms = vec_shared_terms,
                       stringsAsFactors = FALSE)
  return(result)
}
# --- Cluster Matrix ------------------
#' Perform Hierarchical Clustering on a Data Matrix
#'
#' Clusters samples using various distance metrics and clustering approaches to determine optimal cluster numbers.
#'
#' @keywords internal
#' @noRd
.cluster_mat <- function(
    data_matrix,
    dist_method = NULL,
    cluster_method = "ward.D",
    clustering_approach = "gap",
    gap_stat_k_max = 20,
    gap_stat_B = 50,
    density_quantile = 0.90) {

  # Determine appropriate distance metric
  if (is.null(dist_method)) {
    if (all(sapply(data_matrix, is.numeric))) {
      dist_method <- "euclidean"
    } else {
      dist_method <- "gower"
    }
  }

  # Compute distance matrix
  sample_dist <- if (dist_method == "gower") {
    cluster::daisy(data_matrix, metric = "gower")
  } else {
    dist(data_matrix, method = dist_method)
  }

  # Function to determine optimal k based on the selected clustering approach
  determine_k <- function(dist_matrix, cluster_method) {
    if (clustering_approach == "diana") {

      # Determine optimal k using the height difference method
      sample_cluster <- cluster::diana(as.dist(dist_matrix))
      height_diffs <- diff(sample_cluster$height)
      cutoff_index <- which.max(height_diffs)
      return(length(sample_cluster$height) - cutoff_index)

    } else if (clustering_approach == "gap") {

      # Determine optimal k using the gap statistic
      gap_stat <- cluster::clusGap(as.matrix(dist_matrix), FUN = hcut, K.max = gap_stat_k_max, B = gap_stat_B)
      return(cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"]))

    } else if (clustering_approach == "elbow") {

      # Determine optimal k using the elbow method
      wss_plot <- factoextra::fviz_nbclust(as.matrix(dist_matrix), FUN = hcut, method = "wss")

      # Identify the first significant drop and ensure it's a number we can use
      k_optimal <- which.min(diff(diff(wss_plot$data$y))) + 1
      if (is.na(k_optimal) || k_optimal < 2) k_optimal <- 3
      return(k_optimal)

    } else if (clustering_approach == "dynamic") {

    # Determine optimal k using dynamic tree cut
     sample_cluster <- hclust(as.dist(dist_matrix), method = cluster_method)
      cut_clusters <- dynamicTreeCut::cutreeDynamic(
        dendro = sample_cluster,
        distM = as.matrix(as.dist(dist_matrix)),
        deepSplit = 2)
      return(length(unique(cut_clusters)))

    } else if (clustering_approach == "density") {
      # Determine optimal k using density-based clustering
      dclust <- densityClust::densityClust(as.dist(dist_matrix), gaussian = TRUE)
      dclust <- densityClust::findClusters(dclust, rho = quantile(dclust$rho, density_quantile), delta = quantile(dclust$delta, density_quantile))
      return(length(unique(dclust$clusters)))

    } else {
      stop("Invalid clustering approach selected.")
    }
  }

  k_samples <- determine_k(sample_dist, cluster_method)

  # Perform hierarchical clustering only if needed
  sample_cluster <- hclust(as.dist(sample_dist), method = cluster_method)

  # Cut dendrograms using optimal k values
  sample_groups <- cutree(sample_cluster, k = k_samples)

  message("Optimal number of clusters for samples: ", k_samples)

  return(sample_groups)
}
# --- TidyGraph to Data Frame -------------
#' Convert TidyGraph Object to Data Frame
#'
#' Converts a TidyGraph object to a data frame format for easier manipulation and analysis.
#'
#' @keywords internal
#' @noRd
.tidygraph_to_df <- function(tg) {
  require(tidygraph)
  require(dplyr)

  merged_df <- tg |>
    # Grab edge attributes
    tidygraph::activate(edges) |>
    tidygraph::as_tibble() |>

    # Join with node attributes
    dplyr::left_join(
      tg |>
        tidygraph::activate(nodes) |>
        tidygraph::as_tibble() |>
        dplyr::mutate(from = dplyr::row_number()),
      by = "from"
    ) |>

    # Rename columns for clarity
    dplyr::rename(from_name = name) |>

    # Join with node attributes again for the target node
    dplyr::left_join(
      tg |>
        tidygraph::activate(nodes) |>
        tidygraph::as_tibble() |>
        dplyr::mutate(to = dplyr::row_number()),
      by = "to"
    ) |>

    # Rename columns for clarity
    dplyr::rename(to_name = name)

  return(merged_df)
}
# --- Summarize Graph ----------
#' Summarize Graph Attributes
#'
#' Generates a summary of graph attributes including number of nodes, edges, components, diameter, mean distance, and modularity.
#'
#' @keywords internal
#' @noRd
.summarize_graph <- function(g) {
  tibble::tibble(
    "No. of Nodes" = g |> igraph::gorder(),
    "No. of Edges" = g |> igraph::gsize(),
    "No. of Components" = g |> igraph::count_components(),
    "Diameter" = g |> igraph::diameter(),
    "Mean Distance" = g|> igraph::mean_distance(),
    "Modularity" = g |>
      tidygraph::as_tbl_graph() |>
      tidygraph::activate(nodes) |>
      dplyr::mutate(community = tidygraph::group_louvain()) |>
      tidygraph::with_graph(tidygraph::graph_modularity(community))
  )
}
# --- Build Network Plot -------------------
#' Build Network Plot
#'
#' Generates a network plot using ggraph, with options for node color, labels, and facetting.
#'
#' @keywords internal
#' @noRd
.build_ggraph_plot <- function(g,
                               node_color_var,
                               node_colors,
                               label,
                               facet_var,
                               include_stats,
                               fg_text_colour,
                               foreground,
                               alpha,
                               size_lab,
                               color_lab) {
  require(ggplot2)
  require(tidygraph)

  if(!is.null(node_color_var)){
    p <- g |>
      mutate(centrality = tidygraph::centrality_degree()) |>
      ggraph(layout = 'kk') +
      geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
      geom_node_point(aes(size = centrality, color = !!sym(node_color_var))) +
      theme_graph(fg_text_colour = fg_text_colour) +
      labs(size = size_lab,
           color = color_lab)
  } else {
    p <- g |>
      mutate(centrality = tidygraph::centrality_degree()) |>
      ggraph(layout = 'kk') +
      geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
      geom_node_point(aes(size = centrality)) +
      theme_graph(fg_text_colour = fg_text_colour) +
      labs(size = size_lab)
  }

  if(!is.null(node_colors)) {
    p <- p + scale_color_manual(values = node_colors)
  } else {
    n_colors <- g |>
      tidygraph::activate(nodes) |>
      tidygraph::as_tibble() |>
      dplyr::pull(group) |>
      unique() |>
      length()
    p <- p + scale_color_tidy_exp()
  }

  if (!is.null(facet_var)) {
    p <- p +
      ggh4x::facet_grid2(
        as.formula(paste("~", facet_var)),
        scales = "free_x",
        space = "free_x",
        strip = ggh4x::strip_themed(
          background_x = ggh4x::elem_list_rect(
            fill = scales::alpha(foreground, alpha = alpha)
          )
        )
      ) +
      theme(strip.text.x = element_text(face = "bold.italic"))
  }

  if (label) {
    p <- p + geom_node_label(aes(label = label),fontface = "italic", repel = T)
  }

  if (include_stats) {
    summary_info <- .summarize_graph(g)
    label_text <- paste0(
      "Nodes: ", summary_info$`No. of Nodes`,
      " | Edges: ", summary_info$`No. of Edges`,
      " | Components: ", summary_info$`No. of Components`,
      " | Diameter: ", summary_info$Diameter,
      " | Mean Distance: ", round(summary_info$`Mean Distance`, 2),
      " | Modularity: ", round(summary_info$Modularity, 3)
    )

    p <- p +
      labs(caption = label_text) +
      theme(plot.caption = element_text(face = "italic"))
  }

  return(p)
}
# --- Ontology Root --------------
#' Collapse Ontology Root Nodes
#'
#' Collapses a set of ontology terms to either user specified root nodes or root level depending on the level of detail they want
#'
#' @keywords internal
#' @noRd
.collapse_ont_terms <- function(
    ontology_df,
    term_ids,
    root_level = "top",
    assign_label = TRUE
) {
  # Sanitize input: fix character → list columns
  fix_rel_cols <- function(x) {
    if (!is.list(x)) strsplit(x, ";\\s*") else x
  }

  ontology_df <- ontology_df |>
    dplyr::mutate(
      ancestors = fix_rel_cols(ancestors),
      parents = fix_rel_cols(parents),
      children = fix_rel_cols(children)
    )

  get_ancestors <- function(term_id) {
    idx <- match(term_id, ontology_df$id)
    if (!is.na(idx)) ontology_df$ancestors[[idx]] else character(0)
  }

  get_name <- function(term_id) {
    idx <- match(term_id, ontology_df$id)
    if (!is.na(idx)) ontology_df$name[[idx]] else NA_character_
  }

  get_depth <- function(term_id) {
    length(get_ancestors(term_id))
  }

  # Root nodes depend on input
  if (identical(root_level, "top")) {
    root_nodes <- ontology_df$id[
      sapply(ontology_df$parents,
             function(p) length(p) == 0)]

    assign_to_root <- function(term_id) {
      anc <- c(term_id, get_ancestors(term_id))
      matched <- intersect(anc, root_nodes)

      if (length(matched)) matched[1] else NA_character_
    }
  } else if (is.numeric(root_level)) {

    all_depths <- sapply(ontology_df$id, get_depth)
    id_to_depth <- setNames(all_depths, ontology_df$id)

    assign_to_root <- function(term_id) {
      anc <- c(term_id, get_ancestors(term_id))
      if (length(anc) == 0) return(NA_character_)
      # Filter to ancestors that are in the ontology and have depth
      anc_depths <- id_to_depth[anc]
      anc_depths <- anc_depths[!is.na(anc_depths)]
      if (length(anc_depths) == 0) return(NA_character_)

      # Find the deepest matching node ≤ root_level
      eligible <- anc_depths[anc_depths <= root_level]
      if (length(eligible)) {
        return(names(eligible)[which.max(eligible)])  # closest from below
      } else {
        return(NA_character_)
      }
    }
  } else if (is.character(root_level)) {
    root_nodes <- root_level
    assign_to_root <- function(term_id) {
      anc <- c(term_id, get_ancestors(term_id))
      matched <- intersect(anc, root_nodes)
      if (length(matched)) matched[1] else NA_character_
    }
  } else {
    stop("`root_level` must be 'top', a depth integer, or vector of term IDs.")
  }

  assigned_root_ids <- sapply(term_ids, assign_to_root)

  if (assign_label) {
    tibble::tibble(
      term_id = term_ids,
      term_label = sapply(term_ids, get_name),
      root_id = assigned_root_ids,
      root_label = sapply(assigned_root_ids, get_name)
    )
  } else {
    assigned_root_ids
  }
}
