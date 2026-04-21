# Plot Top Features by Factor from Integration Results

Visualizes the top loading features for each factor from multi-omics
integration results (e.g., MOFA, MCIA, DIABLO, RGCCA).

## Usage

``` r
plot_top_factor_features(
  exposomicset,
  feature_col = "feature",
  factors = NULL,
  top_n = 5,
  facet_cols = NULL,
  exp_name_cols = NULL,
  alpha = 0.5
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing integration results in the
  `metadata` slot (must include `integration_results`).

- feature_col:

  A character string indicating the column name to use for y-axis
  feature labels (e.g., `"feature"`, `"gene_symbol"`). This should match
  a column in the output of
  [`pivot_feature()`](https://BioNomad.github.io/tidyexposomics/reference/pivot_feature.md).
  Default is `"feature"`.

- factors:

  Character vector of factors to include (e.g., "Factor1", "Factor2").
  If `NULL`, all factors are plotted.

- top_n:

  Integer specifying the number of top features to show per factor.
  Default is `5`.

- facet_cols:

  Optional color palette for facet strip backgrounds (one per
  `exp_name`), used to distinguish factors.

- exp_name_cols:

  Optional color palette for experiment labels in the plot (`exp_name`),
  passed to `scale_color_manual()`.

- alpha:

  Numeric value between 0 and 1 controlling the transparency of facet
  strip background fill. Default is `0.5`.

## Value

A `ggplot2` object with one facet per factor, showing the top features
and their loadings by experiment.

## Details

This function supports the following integration methods:

- `"MOFA"`: Uses feature weights from MOFA2 (`get_weights()`).

- `"MCIA"`: Uses block loadings from MCIA (`@block_loadings`).

- `"DIABLO"`: Extracts block-specific loadings from `loadings`.

- `"RGCCA"`: Extracts block-specific loadings from `a`.

For each factor, it:

- Selects the top `top_n` features by **absolute loading**.

- Merges with feature metadata using
  [`pivot_feature()`](https://BioNomad.github.io/tidyexposomics/reference/pivot_feature.md).

- Creates a point-range plot showing the loading magnitude.

- Facets each factor with a customizable strip background.

The `feature_col` argument allows you to control which feature-level
metadata column (e.g., gene symbols, metabolite names) is used for
labeling the y-axis.

If palettes are not provided, defaults are chosen using
[`ggpubr::get_palette()`](https://rpkgs.datanovia.com/ggpubr/reference/get_palette.html).

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

mae <- run_multiomics_integration(
    mae,
    method = "DIABLO",
    outcome = "smoker",
    n_factors = 3
)
#> Scaling each assay in MultiAssayExperiment.
#> Running multi-omics integration using DIABLO...
#> Applying DIABLO supervised integration.
#> Design matrix has changed to include Y; each block will be
#>             linked to Y.

# plot top features using default `feature` column
top_feature_p <- mae |>
    plot_top_factor_features()
```
