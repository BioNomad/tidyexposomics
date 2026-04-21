# Plot Circular Network of Exposure Relationships

Generates a circular network plot to visualize relationships between
exposures, either based on correlation ("exposures") or shared features
("degs", "factors").

## Usage

``` r
plot_circos_correlation(
  exposomicset,
  feature_type = c("degs", "omics", "factors", "factor_features", "exposures", "pcs"),
  exposure_cols = NULL,
  corr_threshold = NULL,
  shared_cutoff = 10,
  annotation_colors = NULL,
  low = "#006666",
  mid = "white",
  high = "#8E0152",
  midpoint = NULL
)
```

## Arguments

- exposomicset:

  A MultiAssayExperiment object.

- feature_type:

  One of "exposures", "degs", or "factors".

- exposure_cols:

  Character vector of exposures to include (only for "exposures").

- corr_threshold:

  Minimum \|correlation\| (only for "exposures").

- shared_cutoff:

  Minimum number of shared features (only for "degs" or "factors").
  Default = 10.

- annotation_colors:

  Optional named vector of colors for categories.

- low:

  low value color for edges.

- mid:

  middle value color for edges.

- high:

  high value color for edges.

- midpoint:

  Midpoint for edge color gradient. Defaults to 0 (for correlations) or
  mean shared features.

## Value

A ggplot object (ggraph circular plot).

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

# run correlation analysis
mae <- mae |>
    run_correlation(
        feature_type = "exposures",
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )

# create the circos plot
circos_plot <- mae |>
    plot_circos_correlation(
        feature_type = "exposures"
    )
```
