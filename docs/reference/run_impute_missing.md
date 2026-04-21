# Impute Missing Exposure and Omics Data in a MultiAssayExperiment

Performs missing data imputation on both exposure variables (from
`colData`) and omics datasets (from `experiments`) within a
`MultiAssayExperiment` object.

## Usage

``` r
run_impute_missing(
  exposomicset,
  exposure_impute_method = "median",
  exposure_cols = NULL,
  omics_impute_method = NULL,
  omics_to_impute = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposures and omics data.

- exposure_impute_method:

  Character. Imputation method to use for exposure variables. Defaults
  to `"median"`.

- exposure_cols:

  Character vector. Names of columns in `colData` to impute. If `NULL`,
  all numeric columns are used.

- omics_impute_method:

  Character. Imputation method to use for omics data. Defaults to
  `"knn"`.

- omics_to_impute:

  Character vector. Names of omics datasets to impute. If `NULL`, all
  omics datasets are included.

## Value

A `MultiAssayExperiment` object with imputed exposure and/or omics data.

## Details

For exposures, numeric columns in `colData` are imputed using the
selected method. For omics data, assays are selected and imputed
individually.

Supported imputation methods include:

- `"median"`: Median imputation using
  [`naniar::impute_median_all`](https://naniar.njtierney.com/reference/scoped-impute_median.html)

- `"mean"`: Mean imputation using
  [`naniar::impute_mean_all`](https://naniar.njtierney.com/reference/scoped-impute_mean.html)

- `"knn"`: k-nearest neighbor imputation using
  [`impute::impute.knn`](https://rdrr.io/pkg/impute/man/impute.knn.html)

- `"mice"`: Multiple imputation using chained equations
  ([`mice::mice`](https://amices.org/mice/reference/mice.html))

- `"missforest"`: Random forest-based imputation using
  [`missForest::missForest`](https://rdrr.io/pkg/missForest/man/missForest.html)

- `"lod_sqrt2"`: Substitution of missing values with LOD/sqrt(2), where
  LOD is the smallest non-zero value per variable

## Examples

``` r
# Create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# Introduce some missingness
MultiAssayExperiment::colData(mae)$exposure_pm25[sample(1:20, 5)] <- NA

# Filter features and exposures with high missingness
mae <- run_impute_missing(
    exposomicset = mae,
    exposure_impute_method = "median"
)
#> Imputing exposure data using method: median
```
