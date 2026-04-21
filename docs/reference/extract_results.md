# Extract Results from `MultiAssayExperiment` Metadata

Retrieves a specific analysis result from the metadata slot of a
`MultiAssayExperiment` object.

## Usage

``` r
extract_results(
  exposomicset,
  result = c("codebook", "quality_control", "correlation", "association",
    "differential_analysis", "multiomics_integration", "network", "enrichment")
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object.

- result:

  A character string indicating which result to extract from metadata.
  Must be one of: `"codebook"`, `"quality_control"`, `"correlation"`,
  `"association"`, `"mixture_analysis"`, `"differential_analysis"`,
  `"multiomics_integration"`, `"network"`, `"enrichment"`.

## Value

The corresponding result object stored in `metadata(exposomicset)`, or
`NULL` if not present.

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

# extract results
res <- extract_results(
    exposomicset = mae,
    result = "codebook"
)
```
