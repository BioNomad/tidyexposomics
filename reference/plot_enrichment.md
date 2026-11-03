# Plot Enrichment Results from exposomicset

Visualize enrichment results stored in a `MultiAssayExperiment` object.
Supports dotplots, heatmaps, cnetplots, networks, and multi-panel
summary plots.

## Usage

``` r
plot_enrichment(
  exposomicset,
  feature_type = c("degs", "degs_robust", "omics", "factor_features", "degs_cor",
    "omics_cor", "factor_features_cor"),
  plot_type = c("dotplot", "cnet", "network", "heatmap", "summary"),
  top_n = 5,
  n_per_group = 5,
  add_top_genes = TRUE,
  top_n_genes = 5,
  fill_data = NULL,
  fill_col = "estimate",
  fill_label = NULL,
  overlap_thresh = 0.2,
  node_radius = 0.2,
  pie_colors = NULL,
  label_top_n = NULL,
  label_colour = "black",
  net_facet_by = NULL,
  max_terms = 30,
  node_size = 1,
  term_node_correction = 0.2,
  gene_node_correction = 3,
  go_groups = NULL,
  layout_algo = "fr",
  edge_alpha = 0.3,
  label_size = 3,
  feature_col = "feature"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with enrichment results added via
  [`run_enrichment()`](https://BioNomad.github.io/tidyexposomics/reference/run_enrichment.md).

- feature_type:

  Character; one of `"degs"`, `"degs_robust"`, `"omics"`,
  `"factor_features"`, `"degs_cor"`, `"omics_cor"`, or
  `"factor_features_cor"`. Defines which enrichment results to use.

- plot_type:

  Type of plot to generate. One of `"dotplot"`, `"cnet"`, `"network"`,
  `"heatmap"`, or `"summary"`.

- top_n:

  Integer; number of top `go_group`s to include (used in `"dotplot"`).
  Default is `5`.

- n_per_group:

  Integer; number of terms per group to plot (used in `"dotplot"`).
  Default is `5`.

- add_top_genes:

  Logical; if `TRUE`, appends top shared genes to dotplot facets.
  Default is `TRUE`.

- top_n_genes:

  Integer; number of top genes to show in each group (used in
  `"dotplot"`). Default is `5`.

- fill_data:

  Optional data frame used to color heatmap tiles and cnet gene nodes.
  Needs `exp_name`, `feature`, and a numeric column named by `fill_col`.
  Build it from any source (differential abundance, exposure-omics
  association, or your own scores). If `NULL`, tiles and nodes are grey.

- fill_col:

  Name of the numeric column in `fill_data`. Default is `"estimate"`.

- fill_label:

  Legend title. Defaults to `fill_col`.

- overlap_thresh:

  Numeric; Jaccard threshold for edges in the network plot. Default is
  `0.2`.

- node_radius:

  Numeric; node size in network plot. Default is `0.2`.

- pie_colors:

  Optional named vector of colors for pie charts (network and cnet).

- label_top_n:

  Integer; number of top nodes to label in network. Default is `NULL`.

- label_colour:

  Color of node labels in network. Default is `"black"`.

- net_facet_by:

  Column used to facet the network plot (e.g., `"category"`). Default is
  `NULL`.

- max_terms:

  Integer; max number of terms to include in the cnet plot. Default is
  `30`.

- node_size:

  Numeric; base node size for cnet plot. Default is `1`.

- term_node_correction:

  Scaling factor for term nodes in cnet plot. Default is `0.2`.

- gene_node_correction:

  Scaling factor for gene nodes in cnet plot. Default is `3`.

- go_groups:

  Optional character vector of GO group names to subset enrichment
  results (all plots). Required for `"heatmap"`.

- layout_algo:

  Graph layout algorithm to use in `"network"` and `"cnet"` plots.
  Default is `"fr"`.

- edge_alpha:

  Transparency of network/cnet plot edges. Default is `0.3`.

- label_size:

  Font size for labels in network and cnet plots. Default is `3`.

- feature_col:

  Column in the feature metadata holding the gene identifiers used in
  [`run_enrichment()`](https://BioNomad.github.io/tidyexposomics/reference/run_enrichment.md).
  Must match the `feature_col` used there. Used to map `fill_data` onto
  enriched genes. Default is `"feature"`.

## Value

A `ggplot` or `patchwork` object corresponding to the requested plot
type.

## Details

This function visualizes results from
[`run_enrichment()`](https://BioNomad.github.io/tidyexposomics/reference/run_enrichment.md)
using one of several plot types:

- `"dotplot"`: Enrichment terms grouped by GO group, colored by
  significance.

- `"heatmap"`: Term - gene matrix with optional fill and shared gene
  highlighting.

- `"network"`: Graph of term overlap based on shared genes, faceted by
  metadata if desired.

- `"cnet"`: Gene - term bipartite graph with gene values and term pie
  slices.

- `"summary"`: Multi-panel figure with GO group ridgeplots, gene counts,
  and Venn diagram.

For `"heatmap"` and `"cnet"`, `fill_data` must have one row per
`exp_name` x `feature`. When a gene symbol maps to several features
within an experiment, the value with the largest magnitude is kept.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 30,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# perform differential abundance analysis
mae <- run_differential_abundance(
    exposomicset = mae,
    formula = ~ smoker + sex,
    abundance_col = "counts",
    method = "limma_voom",
    action = "add"
)
#> Running differential abundance testing.
#> Processing assay: mRNA
#> Warning: All samples appear to belong to the same group.
#> Warning: The `.abundance` argument of `test_differential_abundance()` is deprecated as
#> of tidybulk 2.0.0.
#> ℹ Please use the `abundance` argument instead.
#> ℹ The deprecated feature was likely used in the tidyexposomics package.
#>   Please report the issue at
#>   <https://github.com/BioNomad/tidyexposomics/issues>.
#> =====================================
#> tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance
#> or adjust_abundance have been calculated. Therefore, it is essential to add covariates
#> such as batch effects (if applicable) in the formula.
#> =====================================
#> This message is displayed once per session.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Processing assay: proteomics
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Differential abundance testing completed.

# perform enrichment analysis
mae <- run_enrichment(
    exposomicset = mae,
    feature_type = "degs",
    feature_col = "symbol",
    species = "goa_human",
    deg_logfc_threshold = log2(1),
    deg_pval_col = "P.Value",
    deg_pval_threshold = 0.2,
    action = "add"
)
#> 
#> Cannot access https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz. unknown: Unknown error (probably timeout).
#> GO term data failed to download

# create an enrichment plot
enr_plot <- plot_enrichment(
    exposomicset = mae,
    feature_type = "degs",
    plot_type = "network"
)
```
