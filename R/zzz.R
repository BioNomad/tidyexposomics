#' @keywords internal
#' @importFrom stats IQR as.dist as.formula coef cov cutree dist end
#' @importFrom stats glm hclust median na.omit prcomp qchisq qt quantile
#' @importFrom stats reorder rlnorm rnorm runif sd setNames start var cor binomial p.adjust
#' @importFrom utils head tail
#' @importFrom methods as
#' @importFrom dplyr percent_rank row_number
#' @importFrom tibble tibble

utils::globalVariables(c(
    ".exp_name", ".feature", ".sample",
    "abs_loading", "adj.P.Val", "ancestors", "angle", "avg",
    "Category", "category", "celltype", "centrality",
    "centrality_betweenness", "centrality_closeness", "centrality_degree",
    "centrality_eigen", "children", "Cluster", "cluster", "cols", "community",
    "component", "contrast", "correlation",
    "degree", "direction",
    "edge_color", "edge_weight", "edges", "effect", "effect_consistency",
    "end", "estimate", "exp_label", "exp_name", "exp_name_1", "exp_name_2",
    "exp_name_clean", "exp_name_feature", "exp_name_plot",
    "exposome_score_mean", "exposome_score_median", "exposome_score_quantile",
    "exposome_score_sum", "exposome_score_var", "exposure",
    "FDR", "feature", "fill_value", "Freq", "from",
    "gene", "gene_col", "gene_count", "gene_list", "gene_node", "go_group",
    "GO_group", "group",
    "id", "ids", "index",
    "jaccard",
    "label", "loading", "logFC",
    "mean_betweenness", "mean_closeness", "mean_degree", "mean_eigen",
    "mean_log_p", "metric", "missingness", "model",
    "n_above", "n_cat", "n_category", "n_exp", "n_exp_name", "n_exps",
    "n_omics", "n_terms", "n_values", "n_with_sel", "name", "node_index",
    "nodes", "num_shared",
    "omic", "overlap",
    "p_adjust", "p_value", "p.value", "P.Value", "padj", "parents", "PC1",
    "PC2", "pct_miss", "presence_rate",
    "ref_estimate", "ref_se", "reference", "relation",
    "Sample", "sample_group", "scaled_value", "score", "score_tmp", "se",
    "shared", "source_std", "start", "std.error",
    "target", "target_std", "term", "term_name", "term_node", "thresh_met",
    "to", "total", "total_significant", "type",
    "value", "Var1", "var", "var1", "var1_type", "var2", "var2_type",
    "var_exp", "var_label", "var_max", "var_min", "var_pc", "variable",
    "vertex_id",
    "weight",
    "x",
    "y",
    "z_diff",

    # function names
    ":=", "geom_edge_link", "geom_node_arc_bar",
    "geom_node_label", "geom_node_point", "ggraph", "scale_edge_width"
))

#' Internal - onLoad hook to register www assets
#'
#' Registers \code{inst/app/www} as a Shiny resource path \code{"www"} if present.
#'
#' @return Invisibly returns \code{NULL}
#' @keywords internal
#' @importFrom shiny addResourcePath
.onLoad <- function(libname, pkgname) {
    www <- system.file("app", "www", package = pkgname)
    if (dir.exists(www)) shiny::addResourcePath("www", www)
    invisible(NULL)
}
