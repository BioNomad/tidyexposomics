# Plot Network Graph of Features or Exposures

Visualizes network structures created by
[`run_create_network()`](https://BioNomad.github.io/tidyexposomics/reference/run_create_network.md)
from the metadata of a `MultiAssayExperiment` object. Nodes can
represent features (e.g., genes or factors) or exposures, and edges
represent correlations or shared connections.

## Usage

``` r
plot_network(
  exposomicset,
  network = c("degs", "omics", "factors", "factor_features", "exposures",
    "degs_feature_cor", "omics_feature_cor", "factor_features_feature_cor"),
  include_stats = TRUE,
  nodes_to_include = NULL,
  centrality_thresh = NULL,
  top_n_nodes = NULL,
  cor_thresh = NULL,
  label = FALSE,
  label_top_n = 5,
  nodes_to_label = NULL,
  facet_var = NULL,
  foreground = "steelblue",
  fg_text_colour = "grey25",
  node_colors = NULL,
  node_color_var = NULL,
  alpha = 0.5,
  size_lab = "Centrality",
  color_lab = "Group"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing network results in
  metadata.

- network:

  Character string specifying the network type. One of `"degs"`,
  `"omics"`, `"factors"`, `"factor_features"`, `"exposures"`,
  `"degs_feature_cor"`, `"omics_feature_cor"`,
  `"factor_features_feature_cor"`.

- include_stats:

  Logical; if `TRUE`, include edge weights and node centrality metrics
  in the plot aesthetics. Default is `TRUE`.

- nodes_to_include:

  Optional character vector of node names to include (subset of `name`).

- centrality_thresh:

  Optional numeric threshold to filter nodes by centrality degree.

- top_n_nodes:

  Optional integer to keep only the top N nodes by centrality.

- cor_thresh:

  Optional numeric threshold to filter edges by minimum absolute
  correlation.

- label:

  Logical; whether to label nodes. If `TRUE`, top nodes will be labeled.

- label_top_n:

  Integer; number of top-centrality nodes to label if `label = TRUE`.
  Default is `5`.

- nodes_to_label:

  Optional character vector of specific nodes to label.

- facet_var:

  Optional node metadata column to facet the network layout by.

- foreground:

  Color for node outlines and edge borders. Default is `"steelblue"`.

- fg_text_colour:

  Color of node label text. Default is `"grey25"`.

- node_colors:

  Optional named vector of colors for node groups.

- node_color_var:

  Optional node attribute used for node color mapping.

- alpha:

  Alpha transparency for nodes and edges. Default is `0.5`.

- size_lab:

  Legend title for node size (typically centrality). Default is
  `"Centrality"`.

- color_lab:

  Legend title for node color group. Default is `"Group"`.

## Value

A `ggraph` plot object.

## Details

This function retrieves the stored graph object and optionally filters
or labels nodes based on: centrality, correlation, user input, or
group-specific attributes. It supports layout faceting, custom color
mappings, and highlights highly central nodes.

Large graphs (\> 500 nodes) will prompt the user before plotting.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# perform correlation analyses
# correlate with exposures
mae <- mae |>
    run_correlation(
        feature_type = "omics",
        variable_map = mae |>
            pivot_feature() |>
            dplyr::select(
                variable = .feature,
                exp_name = .exp_name
            ),
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )
#> Log2-Transforming each assay in MultiAssayExperiment.

# create the networks
mae <- mae |>
    run_create_network(
        feature_type = "omics",
        action = "add"
    )
#> Creating network from correlation results.
#> Network added to metadata as: network_omics

# plot the network
network_p <- mae |>
    plot_network(
        network = "omics"
    )
#> Extracting graph.
```
