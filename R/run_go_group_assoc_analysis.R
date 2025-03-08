run_go_group_assoc_analysis <- function(expomicset,
                                        geneset, 
                                        outcome, 
                                        mirna_assays = NULL,
                                        covariates = NULL, 
                                        action = "add") {
  library(MultiAssayExperiment)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  
  enrich_res <- MultiAssayExperiment::metadata(expomicset)$functional_enrichment |> 
    pluck(geneset) |> 
    pluck("enrich_res") |> 
    inner_join(
      MultiAssayExperiment::metadata(expomicset)$functional_enrichment |> 
        pluck(geneset) |> 
        pluck("go_groups"), 
      by = "Description") |> 
    mutate(go_group=paste("Group", go_group, sep="_")) |> 
    filter(!grepl("miRNA", exp_name))
  
  if (!is.null(mirna_assays)) {
    enrich_res <- enrich_res |> 
      filter(!exp_name %in% mirna_assays)
  }
  
  pca_res_all <- enrich_res |>
    dplyr::select(exp_name, Cluster, go_group) |> 
    distinct() |> 
    pmap(function(exp_name, Cluster, go_group) {
      
      # Filter enrichment results
      enrich_df_filt <- enrich_res |> 
        filter(exp_name == !!exp_name,
               Cluster == !!Cluster, 
               go_group == !!go_group)
      
      # Get genes per GO group
      genes_per_go_group <- enrich_df_filt |> 
        pull(geneID) |> 
        str_split("/") |> 
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
      
      # If there are less than 10 genes, skip
      if (length(genes_per_go_group) < 10) {
        message(paste("Skipping", exp_name, Cluster, go_group, "- not enough genes"))
        return(NULL)
      }
      
      # Update assay and colData
      exp <- .update_assay_colData(expomicset, exp_name)
      
      # Get assay data
      assay <- assay(exp)
      
      # Filter assay for selected genes
      assay_filt <- assay[rownames(assay) %in% genes_per_go_group, , drop = FALSE]
      
      # Ensure assay_filt has at least 2 features and variance
      if (nrow(assay_filt) < 2 || all(apply(assay_filt, 1, var) == 0)) {
        message(paste("Skipping", exp_name, Cluster, go_group, "- insufficient variance in features"))
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
      names(pca_res_df) <- c(paste("PC", exp_name, Cluster, go_group, sep="/"), "id_to_map")
      
      # Create dataframe for genes and terms
      pca_genes_terms_df <- data.frame(
        term = paste("PC", exp_name, Cluster, go_group, sep="/"),
        top_terms = top_terms,
        genes = ifelse(length(genes_per_go_group) > 0, paste(genes_per_go_group, collapse = ","), "No genes available")
      )
      
      return(list(pca_res_df=pca_res_df,
                  pca_genes_terms_df=pca_genes_terms_df))
    })
  
  pca_res_all <- pca_res_all[unlist(purrr::map(pca_res_all, ~ !is.null(.x)))]
  
  # Filter list to just the pca_res_df
  pca_res_df <- purrr::map(pca_res_all, ~ .x$pca_res_df)
  
  # Filter list to just the pca_genes_terms_df
  pca_genes_terms_df <- purrr::map(pca_res_all, ~ .x$pca_genes_terms_df) |> 
    bind_rows()
  
  model_data <- purrr::reduce(pca_res_df, inner_join, by = "id_to_map") |> 
    inner_join(colData(expomicset) |>
                 as.data.frame() |> 
                 dplyr::select(all_of(c(outcome, covariates))) |> 
                 rownames_to_column("id_to_map"),
               by = "id_to_map") |> 
    column_to_rownames("id_to_map")
  
  pc_cols <- colnames(model_data)[grepl("PC/", colnames(model_data))]
  
  model_data <- model_data |> 
    mutate(across(c(pc_cols, !!sym(outcome)), ~ as.numeric(scale(.x))))
  
  pc_assoc_res <- purrr::map(
    pc_cols,
    function(pc_col) {
      # Create a placeholder column
      model_data$placeholder_col <- model_data[[pc_col]]
      
      model <- glm(
        as.formula(paste(outcome, "~ placeholder_col+", paste(covariates, collapse = "+"))), 
        data = model_data
      )
      
      return(broom::tidy(model) |> 
               mutate(term = ifelse(term == "placeholder_col", pc_col, term)))
    }
  ) |> 
    bind_rows() |> 
    filter(grepl("PC/", term)) |> 
    inner_join(pca_genes_terms_df, by = c("term" = "term")) |> 
    mutate(outcome=outcome)
  
  # Return or store results based on action argument
  if (action == "add") {
    MultiAssayExperiment::metadata(expomicset)$pc_glm_results <- pc_assoc_res
    return(expomicset)
  } else {
    return(pc_assoc_res)
  }
}


x <- expom_enrich |> 
  run_go_group_assoc_analysis(
    geneset = "deg_exp_cor",
    outcome = "fev_height",
    mirna_assays = c("CD4+ T-cell miRNA", "CD16+ Monocyte miRNA"),
    covariates = c("age","sex","race"),
    action = "get")
