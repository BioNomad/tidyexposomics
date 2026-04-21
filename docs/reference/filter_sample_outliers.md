# Filter Sample Outliers

Removes sample outliers from a `MultiAssayExperiment` object based on
PCA analysis.

## Usage

``` r
filter_sample_outliers(exposomicset, outliers = NULL)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics and exposure data.

- outliers:

  An optional character vector specifying sample names to be removed. If
  `NULL`, the function uses outliers identified in
  `metadata(exposomicset)$pca$outliers`. Default is `NULL`.

## Value

A `MultiAssayExperiment` object with the specified outliers removed.

## Details

The function checks for the presence of PCA results in
`metadata(exposomicset)`. If `outliers` is not provided, it retrieves
precomputed outliers from `metadata(exposomicset)$pca$outliers`. The
identified samples are removed from the dataset.

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

# run PCA
mae <- mae |>
    run_pca()
#> Identifying common samples.
#> Subsetting exposure data.
#> Subsetting omics data.
#> Performing PCA on Feature Space.
#> Performing PCA on Sample Space.
#> No outliers detected.

# filter outliers if present
mae <- mae |>
    filter_sample_outliers()
#> Removing outliers: 
```
