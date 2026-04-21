# Plot Correlation Summary from Exposure-Feature Correlations

Generates a bar plot summary of exposure-feature correlations using
customizable modes.

## Usage

``` r
plot_correlation_summary(
  exposomicset,
  feature_type = c("degs", "omics", "factors", "factor_features", "exposures", "pcs"),
  mode = c("top_exposures", "top_features", "exposure_category", "assay", "summary"),
  top_n = 15
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with correlation results in metadata.

- feature_type:

  One of `"degs"`, `"factors"`, `"omics"`, or `"exposures"`.

- mode:

  One of:

  - `"top_exposures"`: Top exposures by assay

  - `"top_features"`: Top features by exposure category

  - `"exposure_category"`: Total associations by exposure category

  - `"assay"`: Total associations by omics assay

  - `"summary"`: Patchwork layout combining all

- top_n:

  Number of top exposures or features to display (for top modes).
  Default is `15`.

## Value

A `ggplot2` object or a `patchwork` object (if `mode = "summary"`).

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

# create the correlation summary plot
cor_summary_plot <- mae |>
    plot_correlation_summary(
        feature_type = "omics",
        mode = "summary"
    )
```
