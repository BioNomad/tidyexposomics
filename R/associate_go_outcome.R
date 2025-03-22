#' Associate Gene Ontology (GO) Enrichment with an Outcome
#'
#' Performs principal component analysis (PCA) on genes within GO enrichment results
#' and tests for associations with an outcome using generalized linear models (GLMs).
#' Optionally stores results in `expomicset` metadata.
#'
#' @param expomicset A `MultiAssayExperiment` object containing functional enrichment data.
#' @param geneset A character string specifying the name of the gene set to analyze.
#' @param outcome A character string specifying the outcome variable name.
#' @param mirna_assays An optional character vector of miRNA assays to exclude. Default is `NULL`.
#' @param covariates An optional character vector of covariate names. Default is `NULL`.
#' @param action A character string indicating whether to return results (`"get"`) or add them to metadata (`"add"`). Default is `"add"`.
#'
#' @details
#' The function retrieves functional enrichment results from `metadata(expomicset)`, filters based on the selected gene set,
#' and extracts genes for PCA. The first principal component (PC1) of each GO enrichment group is computed and tested
#' for association with the outcome using a GLM. The results include model statistics and gene annotations.
#'
#' @return If `action = "add"`, returns the modified `expomicset` with PCA-GLM results stored in metadata.
#' If `action = "get"`, returns a data frame containing:
#' \item{term}{The PC term tested in the GLM.}
#' \item{estimate}{The regression coefficient for the association.}
#' \item{std.error}{The standard error of the estimate.}
#' \item{statistic}{The test statistic from the GLM.}
#' \item{p.value}{The p-value of the association test.}
#' \item{top_terms}{The top GO terms associated with the PC.}
#' \item{genes}{The genes contributing to the PC.}
#' \item{outcome}{The outcome variable tested.}
#'
#' @examples
#' \dontrun{
#' results <- associate_go_outcome(
#'   expomicset = expom,
#'   geneset = "deg_exp_cor",
#'   outcome = "fev_height",
#'   mirna_assays = c("CD4+ T-cell miRNA"),
#'   covariates = c("age", "sex"),
#'   action = "get"
#' )
#' }
#'
#' @export
associate_go_outcome <- function(expomicset,
                                 geneset,
                                 outcome,
                                 feature_col = NULL,
                                 mirna_assays = NULL,
                                 covariates = NULL,
                                 min_genes = 10,
                                 family = "gaussian",
                                 action = "add") {
  # Check if expomicset is a MultiAssayExperiment
  enrich_res <- MultiAssayExperiment::metadata(expomicset)$functional_enrichment |>
    purrr::pluck(geneset)

  # Check if mirna_assays is NULL
  if (!is.null(mirna_assays)) {
    enrich_res <- enrich_res |>
      dplyr::filter(!exp_name %in% mirna_assays)
  }

  pca_res_all <- enrich_res |>
    dplyr::select(exp_name, Cluster, go_group) |>
    dplyr::distinct() |>
    purrr::pmap(function(exp_name, Cluster, go_group) {

      # Filter enrichment results
      enrich_df_filt <- enrich_res |>
        dplyr::filter(exp_name == !!exp_name,
                      Cluster == !!Cluster,
                      go_group == !!go_group)

      # Get genes per GO group
      genes_per_go_group <- enrich_df_filt |>
        dplyr::pull(geneID) |>
        stringr::str_split("/") |>
        unlist() |>
        unique()

      # Get top 3 GO term descriptions
      top_terms <- enrich_df_filt |>
        dplyr::pull(Description) |>
        unique() |>
        head(3) |>
        paste(collapse = "; ")

      # If no top terms are found, assign a placeholder
      if (is.na(top_terms) || top_terms == "") {
        top_terms <- "No GO terms available"
      }

      # If there are less than min_genes genes, skip
      if (length(genes_per_go_group) < min_genes) {
        message(paste("Skipping", exp_name, Cluster, go_group, "- not enough genes"))
        return(NULL)
      }

      # Update assay and colData
      exp <- .update_assay_colData(expomicset, exp_name)

      # If feature_col is not null then map the features to the genes
      if (!is.null(feature_col)) {
        genes_per_go_group <- exp |>
          tidybulk::pivot_transcript() |>
          filter(!!sym(feature_col) %in% genes_per_go_group) |>
          pull(.feature)
      }

      # Get assay data
      assay <- SummarizedExperiment::assay(exp)

      # Filter assay for selected genes
      assay_filt <- assay[rownames(assay) %in% genes_per_go_group, , drop = FALSE]

      # Ensure assay_filt has at least 2 features and variance
      if (nrow(assay_filt) < 2 || all(apply(assay_filt, 1, var) == 0)) {
        message(paste("Skipping",
                      exp_name,
                      Cluster,
                      go_group,
                      "- insufficient variance in features"))
        return(NULL)
      }

      # Transpose to make samples rows
      assay_filt <- t(assay_filt)

      # Apply log transformation to stabilize variance
      assay_filt <- log2(assay_filt + 1)

      # Perform PCA safely
      pca_res <- prcomp(assay_filt, scale. = TRUE)

      # Create output dataframe
      pca_res_df <- data.frame(
        PC_exp_mat = pca_res$x[,1],
        id_to_map = rownames(pca_res$x)
      )

      # Rename PC1 column to indicate the source
      names(pca_res_df) <- c(paste(
        "PC",
        exp_name,
        Cluster,
        go_group,
        sep="/"),
        "id_to_map")

      # Create dataframe for genes and terms
      pca_genes_terms_df <- data.frame(
        term = paste("PC", exp_name, Cluster, go_group, sep="/"),
        top_terms = top_terms,
        genes = ifelse(length(genes_per_go_group) > 0,
                       paste(genes_per_go_group, collapse = ","),
                       "No genes available")
      )

      return(list(pca_res_df=pca_res_df,
                  pca_genes_terms_df=pca_genes_terms_df))
    })

  pca_res_all <- pca_res_all[unlist(purrr::map(pca_res_all, ~ !is.null(.x)))]

  # Filter list to just the pca_res_df
  pca_res_df <- purrr::map(pca_res_all, ~ .x$pca_res_df)

  # Filter list to just the pca_genes_terms_df
  pca_genes_terms_df <- purrr::map(pca_res_all, ~ .x$pca_genes_terms_df) |>
    dplyr::bind_rows()

  model_data <- purrr::reduce(pca_res_df,
                              inner_join,
                              by = "id_to_map") |>
    dplyr::inner_join(MultiAssayExperiment::colData(expomicset) |>
                        as.data.frame() |>
                        dplyr::select(all_of(c(outcome, covariates))) |>
                        tibble::rownames_to_column("id_to_map"),
                      by = "id_to_map") |>
    tibble::column_to_rownames("id_to_map")

  pc_cols <- colnames(model_data)[grepl("PC/", colnames(model_data))]

  if(family=="gaussian"){
    model_data <- model_data |>
      dplyr::mutate(
        dplyr::across(c(pc_cols, !!sym(outcome)),
                      ~ as.numeric(scale(.x))))
  }else if(family=="binomial"){
    model_data <- model_data |>
      dplyr::mutate(
        dplyr::across(c(pc_cols),
                      ~ as.numeric(scale(.x))))
  }

  pc_assoc_res <- purrr::map(
    pc_cols,
    function(pc_col) {
      # Create a placeholder column
      model_data$placeholder_col <- model_data[[pc_col]]

      model <- glm(
        as.formula(paste(outcome, "~ placeholder_col+", paste(covariates, collapse = "+"))),
        data = model_data,
        family = family
      )

      return(broom::tidy(model) |>
               dplyr::mutate(
                 term = ifelse(term == "placeholder_col",
                               pc_col,
                               term)))
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::filter(grepl("PC/", term)) |>
    dplyr::inner_join(pca_genes_terms_df, by = c("term" = "term")) |>
    dplyr::mutate(outcome=outcome)

  # Return or store results based on action argument
  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$pc_glm_results <- pc_assoc_res
    return(expomicset)
  } else {
    return(pc_assoc_res)
  }
}
