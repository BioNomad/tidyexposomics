#' Create and store a network from correlation results in a MultiAssayExperiment
#'
#' Constructs an undirected network based on correlation results between exposures
#' and omics data, and optionally stores it in the metadata of a `MultiAssayExperiment` object.
#' Nodes are derived from both variables in the correlation results, and grouped by data type
#' (e.g., "Gene Expression", "Metabolomics", etc.).
#'
#' @param expomicset A `MultiAssayExperiment` object that contains correlation results
#'   in its metadata.
#' @param cor_results Character string indicating which correlation result to use.
#'   Options are `"deg_exp_cor"` (default), `"factor_exp_cor"`, or `"omic_exp_cor"`.
#' @param action Character string indicating whether to `"add"` the resulting graph to
#'   the object's metadata or to `"get"` it directly. Default is `"add"`.
#'
#' @return If `action = "add"`, returns the modified `MultiAssayExperiment` with the
#'   network added to its metadata. If `action = "get"`, returns a list with the graph
#'   object and a summary of the graph structure.
#'
#' @details
#' The function identifies the appropriate correlation results table based on the
#' `cor_results` argument. It then constructs an `igraph` object from the edges
#' (correlations) and derives node metadata from variable types or categories.
#'
#' The resulting graph is undirected and includes group labels for visualization.
#'
#' @import igraph tidygraph MultiAssayExperiment dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a network and add it to the MultiAssayExperiment metadata
#' mae <- run_create_network(mae, cor_results = "omic_exp_cor", action = "add")
#'
#' # Retrieve the network object without adding it
#' net_obj <- run_create_network(mae, cor_results = "deg_exp_cor", action = "get")
#' igraph::plot(net_obj$graph)
#' }

correlate_exposoures_omics <- function(expomicset,
                                       top_n_omics = NULL,
                                       correlation_method = "spearman",
                                       correlation_cutoff = 0.3,
                                       cor_pval_column = "p.value",
                                       pval_cutoff = 0.05,
                                       action = "add"){
  # Extract and preprocess colData
  message("Extracting exposure data...")
  data <- MultiAssayExperiment::colData(expomicset)  |>
    as.data.frame() |>
    dplyr::select_if(is.numeric)

  # Remove Principal Components if in the meta data
  if(any(grepl("^PC", colnames(data)))){
    data <- data |>
      dplyr::select(-dplyr::matches("^PC"))
  }

  # Log2-transform omics data
  log_2_omics <- .log2_multiassay(expomicset)

  # Check if top_n_omics is provided
  if(!is.null(top_n_omics)){
    selected_features <- .top_var_multiassay(
      expomicset = log_2_omics,
      n = top_n_omics
    )
  } else{
    # Select all features if top_n_omics is not provided
    selected_features <- lapply(
      MultiAssayExperiment::experiments(log_2_omics),
      function(x) {
        rownames(x)
      }
    )
  }


  # Extract log2-transformed omics data
  omics_df <- lapply(names(MultiAssayExperiment::experiments(log_2_omics)),function(name){
    # Extract omics data
    se=MultiAssayExperiment::experiments(log_2_omics)[[name]]

    # Filter omics data based on selected features
    se=se[rownames(se) %in% selected_features[[name]],]

    # Convert assay data to data frame
    df=se |>
      SummarizedExperiment::assay() |>
      t() |>
      as.data.frame()
    names(df)=
      paste0(name,"_",names(df))
    names(df)=gsub(" |-","_",names(df))

    # Create a unique identifier for each row
    df=df |>
      rownames_to_column("id")
  }) |>
    (\(lst){
      # Name list
      names(lst) <- names(MultiAssayExperiment::experiments(log_2_omics))
      lst
    })() |>
    # Combine all omics data into a single data frame
    purrr::reduce(full_join, by = "id")

  # Merge omics data with exposure data
  merged_data <- data |>
    dplyr::mutate(id = rownames(data)) |>
    dplyr::left_join(omics_df, by = "id")

  # Ensure there are no NA values
  merged_data <- na.omit(merged_data)

  if(ncol(data)>5000){
    message("Warning: The number of features is greater than 5000. This may take a while.")
    message("Consider using the top_n_omics parameter to reduce the number of features.")
  }

  # Perform correlation analysis
  message("Running Correlation analysis...")
  correlation_matrix <- merged_data |>
    mutate_all(as.numeric) |>
    as.matrix() |>
    Hmisc::rcorr(type = correlation_method)

  # Convert correlation and p-values to tidy format
  correlation_df <- correlation_matrix$r |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    reshape2::melt(id.vars = "var1") |>
    `colnames<-`(c("var1", "var2", "correlation"))

  pvalue_df <- correlation_matrix$P |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    reshape2::melt(id.vars = "var1") |>
    `colnames<-`(c("var1", "var2", "p.value"))

  # Merge correlation results
  correlation_results <- correlation_df |>
    dplyr::inner_join(pvalue_df, by = c("var1", "var2")) |>
    dplyr::filter(abs(correlation) > correlation_cutoff) |>
    dplyr::mutate(FDR = p.adjust(p.value, method = "fdr")) |>
    dplyr::filter(!!sym(cor_pval_column) < pval_cutoff) |>
    dplyr::arrange(desc(abs(correlation))) |>
    dplyr::mutate(
      var1 = as.character(var1),
      var2 = as.character(var2),
      var_min = pmin(var1, var2),
      var_max = pmax(var1, var2)
    ) |>
    dplyr::select(var_min, var_max, correlation, p.value, FDR) |>
    dplyr::distinct() |>
    dplyr::rename(var1 = var_min, var2 = var_max) |>
    dplyr::filter(var1 != var2)

  # Add annotation columns
  classify_variable <- function(varname) {
    des <- MultiAssayExperiment::metadata(expomicset)$var_info

    if (grepl("^Gene_Expression_", varname)) return("Gene Expression")
    if (grepl("^Proteomics_", varname)) return("Proteomics")
    if (grepl("^Metabolomics_", varname)) return("Metabolomics")
    if (varname %in% des$variable) {
      matched_cat <- as.character(des$category[des$variable == varname])
      return(matched_cat)
    }
    return("Unknown")
  }

  correlation_results <- correlation_results |>
    dplyr::mutate(
      var1_type = purrr::map_chr(var1, classify_variable),
      var2_type = purrr::map_chr(var2, classify_variable)
    )

  message("Correlation analysis completed.")

  if(action=="add"){
    # Save results in metadata
    MultiAssayExperiment::metadata(expomicset)$omics_exposure_correlation <- correlation_results

    return(expomicset)
  }else if (action=="get"){
    return(correlation_results)
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }

}
