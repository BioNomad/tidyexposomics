# List Stored Association Models

Summarise every association model stored in a `MultiAssayExperiment`,
across all sources.

## Usage

``` r
association_models(exposomicset)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` with results from
  [`run_association()`](https://BioNomad.github.io/tidyexposomics/reference/run_association.md).

## Value

A tibble with one row per stored model: `source`, `model`, `formula`,
`subset`, `stratum`, `family`, `n_samples`, and `n_features`.

## Examples

``` r
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
) |>
    run_association(
        formula = smoker ~ age,
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial"
    )
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.
#> Running GLMs.

association_models(mae)
#> # A tibble: 1 × 8
#>   source    model   formula      subset stratum family   n_samples n_features
#>   <chr>     <chr>   <chr>        <chr>  <chr>   <chr>        <int>      <int>
#> 1 exposures default smoker ~ age NA     NA      binomial        10          2
```
