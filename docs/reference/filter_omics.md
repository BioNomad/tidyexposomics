# Filter low-quality features in omics assays

This function applies variance- or expression-based filtering across one
or more assays within a `MultiAssayExperiment` object. It is useful for
removing low-quality or uninformative features before downstream
analysis.

## Usage

``` r
filter_omics(
  exposomicset,
  method = c("variance", "expression"),
  assays = NULL,
  assay_name = 1,
  min_var = 1e-05,
  min_value = 5,
  min_prop = 0.7,
  verbose = TRUE
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics assays.

- method:

  Filtering method: either `"variance"` or `"expression"`.

- assays:

  Character vector of assay names to filter. If `NULL`, all assays are
  filtered.

- assay_name:

  Name or index of the assay within each `SummarizedExperiment` to use.

- min_var:

  Minimum variance threshold (used if `method = "variance"`).

- min_value:

  Minimum expression value (used if `method = "expression"`).

- min_prop:

  Minimum proportion of samples exceeding `min_value` (used if
  `method = "expression"`).

- verbose:

  Whether to print messages for each assay being filtered.

## Value

A filtered `MultiAssayExperiment` object with updated assays and step
record.

## Examples

``` r
# Filter the proteomics assay by variance
filtered_mae <- filter_omics(
    exposomicset = make_example_data(return_mae = TRUE),
    method = c("variance"),
    assays = "proteomics",
    assay_name = 1,
    min_var = 0.01,
    verbose = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.
#> Filtering assay: proteomics
#> Filtered 0 of 80 features from 'proteomics' using method 'variance'
```
