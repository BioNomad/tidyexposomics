# Extract Feature Metadata from a MultiAssayExperiment

Extracts feature-level metadata across all assays in a
`MultiAssayExperiment` and returns a combined tibble.

## Usage

``` r
pivot_feature(exposomicset)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object.

## Value

A tibble with feature metadata from all assays, with an added
`.exp_name` column.

## Details

This function:

- Iterates over all assays in the `MultiAssayExperiment`.

- Updates each assay's sample metadata (`colData`) using
  `.update_assay_colData()`.

- Extracts feature-level metadata using `pivot_transcript()` from
  `tidybulk`.

- Combines results across assays into a single tibble, adding a
  `.exp_name` column.

## Examples

``` r
#' # create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# pivot experiment
feature_data <- mae |>
    pivot_feature()
```
