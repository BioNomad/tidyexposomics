# Generate Example Data for Testing

This helper function generates a reproducible dummy dataset containing
exposures, mRNA data, and proteomics data. It can optionally return the
data as a `MultiAssayExperiment` using
[`create_exposomicset`](https://BioNomad.github.io/tidyexposomics/reference/create_exposomicset.md).

## Usage

``` r
make_example_data(
  n_samples = 12,
  n_proteins = 80,
  use_batch = FALSE,
  return_mae = FALSE
)
```

## Arguments

- n_samples:

  Integer. Number of samples to simulate (default: 12).

- n_proteins:

  Integer. Number of proteins to simulate (default: 80).

- use_batch:

  Logical. If `TRUE`, include a "batch" variable in the exposure data
  (default: FALSE).

- return_mae:

  Logical. If `TRUE`, return a `MultiAssayExperiment` created using
  [`create_exposomicset()`](https://BioNomad.github.io/tidyexposomics/reference/create_exposomicset.md)
  (default: FALSE).

## Value

Either:

- A named list containing `codebook`, `exposure`, `omics`, and
  `row_data`, if `return_mae = FALSE`.

- A `MultiAssayExperiment`, if `return_mae = TRUE`.

## Examples

``` r
# Return as a list
dummy <- make_example_data()

# Return as a MultiAssayExperiment
mae <- make_example_data(return_mae = TRUE)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.
```
