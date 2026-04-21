# Extract Sample Metadata from MultiAssayExperiment or SummarizedExperiment

Extracts and formats sample-level metadata (`colData`) from a
`MultiAssayExperiment` or `SummarizedExperiment` object.

## Usage

``` r
pivot_sample(x, ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` or `SummarizedExperiment` object.

- ...:

  Additional arguments passed to `pivot_sample()` from `tidybulk` for
  `SummarizedExperiment` objects.

## Value

A tibble containing sample metadata with an added `.sample` column.

## Details

This function:

- Extracts **sample metadata** from `MultiAssayExperiment` using
  `colData()`, converting it to a tibble.

- Calls `pivot_sample()` from `tidybulk` when applied to a
  `SummarizedExperiment` object.

- **Error Handling**: Returns an error if `x` is not a
  `MultiAssayExperiment` or `SummarizedExperiment`.

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

sample_data <- mae |>
    pivot_sample()
```
