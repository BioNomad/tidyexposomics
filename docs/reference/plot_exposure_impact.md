# Plot Exposure Impact on Network Centrality Metrics

Visualizes the impact of exposures on network centrality measures of
associated features (e.g., genes or latent factors) as a heatmap. Each
exposure is scored by four centrality metrics, scaled within metric, and
grouped by exposure category.

## Usage

``` r
plot_exposure_impact(
  exposomicset,
  feature_type = c("degs", "omics", "factors"),
  min_per_group = 5,
  facet_cols = NULL,
  bar_cols = NULL,
  alpha = 0.3,
  ncol = 2,
  nrow = 1,
  heights = c(1, 1),
  widths = c(2, 1)
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with results from
  [`run_exposure_impact()`](https://BioNomad.github.io/tidyexposomics/reference/run_exposure_impact.md).

- feature_type:

  Character string specifying the feature type. One of `"degs"`,
  `"omics"`, or `"factors"`.

- min_per_group:

  Minimum number of features per exposure for inclusion (not currently
  used). Default is `5`.

- facet_cols:

  Optional named vector of colors for exposure categories.

- bar_cols:

  Optional vector of colors for bar plots (if enabled).

- alpha:

  Transparency level for category strips (if enabled). Default is `0.3`.

- ncol, nrow:

  Layout for optional patchwork combination (currently unused). Default:
  `ncol = 2`, `nrow = 1`.

- heights:

  Relative heights and widths for combined plots (currently unused).
  Default: `c(1,1)`.

- widths:

  Relative widths for combined plots (currently unused). Default:
  `c(2,1)`.

## Value

A `ggplot`/`patchwork` object showing a heatmap of scaled network
centrality scores per exposure, annotated by category.

## Details

This function uses the output of
[`run_exposure_impact()`](https://BioNomad.github.io/tidyexposomics/reference/run_exposure_impact.md)
to calculate and visualize the mean centrality values for each exposure
across its associated features. The following network centrality metrics
are shown:

- Degree centrality

- Eigenvector centrality

- Closeness centrality

- Betweenness centrality

All values are scaled within metric across exposures. A side bar
indicates the category of each exposure.

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
    ) |>
    run_correlation(
        feature_type = "omics",
        variable_map = mae |>
            pivot_feature() |>
            dplyr::select(
                variable = .feature,
                exp_name = .exp_name
            ),
        feature_cors = TRUE,
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )
#> Log2-Transforming each assay in MultiAssayExperiment.
#> Log2-Transforming each assay in MultiAssayExperiment.

# create the networks
mae <- mae |>
    run_create_network(
        feature_type = "omics_feature_cor",
        action = "add"
    ) |>
    run_create_network(
        feature_type = "omics",
        action = "add"
    )
#> Creating network from correlation results.
#> Network added to metadata as: network_omics_feature_cor
#> Creating network from correlation results.
#> Network added to metadata as: network_omics

# perform impact analysis
mae <- mae |>
    run_exposure_impact(
        feature_type = "omics"
    )

# create the exposure impact plot
exposure_impact_p <- mae |>
    plot_exposure_impact(
        feature_type = "omics"
    )
```
