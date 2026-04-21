# Perform Principal Component Analysis (PCA)

Runs PCA on the feature and sample spaces of a `MultiAssayExperiment`
object, identifying outliers based on Mahalanobis distance.

## Usage

``` r
run_pca(
  exposomicset,
  log_trans_exp = FALSE,
  log_trans_omics = TRUE,
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics and exposure data.

- log_trans_exp:

  A boolean value specifying whether to log2 transform the exposure data

- log_trans_omics:

  a boolean value specifying whether to log2 transform the omics data

- action:

  A character string specifying whether to store (`"add"`) or return
  (`"get"`) the results. Default is `"add"`.

## Value

If `action = "add"`, the function returns the input
`MultiAssayExperiment` with:

- **PC scores** added as columns in `colData(exposomicset)`, and

- **PCA objects** stored under
  `metadata(exposomicset)$quality_control$pca`.

If `action = "get"`, the function returns a list containing:

- pca_df:

  A tibble of the transformed input data.

- pca_feature:

  A `prcomp` object containing PCA results for features.

- pca_sample:

  A `prcomp` object containing PCA results for samples.

- outliers:

  A character vector of detected sample outliers.

## Details

This function:

- Identifies **common samples** across all assays and exposure data.

- Performs **PCA on features** (transformed and standardized).

- Performs **PCA on samples** and computes Mahalanobis distance to
  detect outliers.

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
```
