# Plot Correlation Tilemap

Visualizes a correlation matrix as a heatmap tile plot using correlation
results stored in the metadata of a `MultiAssayExperiment` object. When
`feature_type = "pcs"`, the function forces PCs to appear on the x-axis
and exposures on the y-axis, and it adds a barplot showing how many PCs
are significantly associated with each exposure. It also suppresses
nonsignificant tiles based on a specified p-value threshold.

## Usage

``` r
plot_correlation_tile(
  exposomicset,
  feature_type = c("pcs", "degs", "omics", "factors", "factor_features", "exposures"),
  pval_cutoff = 0.05,
  na_color = "grey100",
  fill_limits = c(-1, 1),
  midpoint = 0
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing correlation results in
  metadata.

- feature_type:

  Type of correlation results to plot. One of `"pcs"`, `"degs"`,
  `"omics"`, `"factors"`, `"factor_features"`, or `"exposures"`. Must
  match the key used in
  `metadata(exposomicset)$correlation[[feature_type]]`.

- pval_cutoff:

  Numeric p-value cutoff below which correlations are displayed.
  Nonsignificant tiles are rendered in the `na_color`. Default is
  `0.05`.

- na_color:

  Color used to represent nonsignificant or missing correlations.
  Default is `"grey100"`.

- fill_limits:

  Numeric vector of length 2 defining the scale limits for correlation
  values. Default is `c(-1, 1)`.

- midpoint:

  Numeric value for centering the fill gradient. Default is `0`.

## Value

A `ggplot2` tile plot (or a combined tile + barplot if
`feature_type = "pcs"`).

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

# run pca
mae <- mae |>
    run_pca()
#> Identifying common samples.
#> Subsetting exposure data.
#> Subsetting omics data.
#> Performing PCA on Feature Space.
#> Performing PCA on Sample Space.
#> No outliers detected.

# correlate with exposures
mae <- mae |>
    run_correlation(
        feature_type = "pcs",
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )

# make the correlation tile plot
cor_tile_p <- mae |>
    plot_correlation_tile(
        feature_type = "pcs"
    )
```
