# Extract Merged Omics and Exposure Data Frame

This function extracts and merges exposure variables from `colData` with
selected features from omics datasets in a `MultiAssayExperiment`
object. Optionally applies log2 transformation to omics data and
restricts features based on a variable map.

## Usage

``` r
extract_omics_exposure_df(exposomicset, variable_map = NULL, log2_trans = TRUE)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics and exposure data.

- variable_map:

  A data frame with columns `"variable"` and `"exp_name"`, indicating
  which variables belong to each omics or exposure domain.

- log2_trans:

  Logical; whether to log2-transform omics data. Default is `TRUE`.

## Value

A data frame where rows correspond to samples, and columns contain
exposure variables and log2-transformed omics features. Columns from
different omics types are disambiguated using prefixes.

## Details

If `variable_map` is provided, it is used to select variables from both
exposures and omics. If not provided, all numeric `colData` variables
are used as exposures (excluding variables matching `^PC`), and all
omics features are included.

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
# export the omics exposure df
merged_df <- extract_omics_exposure_df(
    mae,
    log2_trans = TRUE
)
#> Log2-Transforming each assay in MultiAssayExperiment.
```
