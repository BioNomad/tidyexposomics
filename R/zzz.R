#' @keywords internal
#' @importFrom stats IQR as.dist as.formula coef cov cutree dist end
#' @importFrom stats glm hclust median na.omit prcomp qchisq qt quantile
#' @importFrom stats reorder rlnorm rnorm runif sd setNames start var
#' @importFrom utils head tail
#' @importFrom methods as

utils::globalVariables(c(
    ".exp_name", ".feature", ".sample", ":=", "binomial", "contrast",
    "Category", "Cluster", "FDR", "Freq", "GO_group", "PC1", "PC2", "Sample",
    "Var1", "abs_loading", "ancestors", "angle", "avg", "category",
    "centrality", "centrality_betweenness", "centrality_closeness",
    "centrality_degree", "centrality_eigen", "children", "cluster", "cols",
    "community", "component", "correlation", "degree", "direction", "edges",
    "edge_weight", "effect_consistency", "end", "estimate", "exp_label",
    "exp_name", "exp_name_1", "exp_name_2", "exp_name_clean",
    "exp_name_feature", "exp_name_plot", "exposure", "exposome_score_mean",
    "exposome_score_median", "exposome_score_quantile", "exposome_score_sum",
    "exposome_score_var", "feature", "from", "gene_col", "gene_count",
    "gene_list", "gene_node", "geom_edge_link", "geom_node_arc_bar",
    "geom_node_label", "geom_node_point", "ggraph", "go_group", "group",
    "id", "ids", "index", "jaccard", "label", "loading", "mean_betweenness",
    "mean_closeness", "mean_degree", "mean_eigen", "mean_log_p", "metric",
    "missingness", "n_above", "n_cat", "n_category", "n_exp", "n_exp_name",
    "n_exps", "n_omics", "n_terms", "n_values", "n_with_sel", "name",
    "node_index", "nodes", "num_shared", "omic", "overlap", "p.adjust",
    "p.value", "p_adjust", "p_value", "padj", "parents", "pct_miss",
    "presence_rate", "reorder", "sample_group", "scale_edge_width",
    "scaled_value", "score", "score_tmp", "se", "setNames", "shared",
    "source_std", "start", "std.error", "target", "target_std", "term",
    "term_name", "term_node", "thresh_met", "to", "total",
    "total_significant", "type", "value", "var", "var1", "var1_type", "var2",
    "var2_type", "var_exp", "var_label", "var_max", "var_min", "var_pc",
    "variable", "weight", "x", "y"
))

#' Internal - onLoad hook to register www assets
#'
#' Registers \code{inst/app/www} as a Shiny resource path \code{"www"} if present.
#'
#' @keywords internal
#' @importFrom shiny addResourcePath
.onLoad <- function(libname, pkgname) {
  www <- system.file("app", "www", package = pkgname)
  if (dir.exists(www)) shiny::addResourcePath("www", www)
}
