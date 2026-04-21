# Plot Summary of Factor Contributions from Multi-Omics Integration

Generates a summary plot of factor contributions from multi-omics
integration results stored in a `MultiAssayExperiment` object.

## Usage

``` r
plot_factor_summary(
  exposomicset,
  low = "#006666",
  mid = "white",
  high = "#8E0152",
  midpoint = 0.5
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing integration results in
  `metadata(exposomicset)$multiomics_integration$integration_results`.

- low:

  Color for low values in the fill gradient. Default is `"#006666"`.

- mid:

  Color for midpoint in the fill gradient. Default is `"white"`.

- high:

  Color for high values in the fill gradient. Default is `"#8E0152"`.

- midpoint:

  Midpoint value for the gradient color scale. Default is `0.5`.

## Value

A `ggplot` object showing factor contributions based on the integration
method.

## Details

This function visualizes factor contributions based on the integration
method:

- **MOFA**: Variance explained per factor and view.

- **MCIA**: Block score weights per omic.

- **DIABLO**: Mean absolute sample score per omic and factor (from
  block-specific variates).

- **RGCCA**: Mean absolute sample score per omic and factor (from
  aligned block scores).

The color gradient can be customized using the `low`, `mid`, `high`, and
`midpoint` parameters.

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

factor_sum_plot <- mae |>
    plot_factor_summary()
```
