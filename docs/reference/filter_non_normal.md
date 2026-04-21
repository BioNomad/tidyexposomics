# Filter Non-Normal Exposure Variables

Removes exposure variables that deviate significantly from a normal
distribution based on normality test results stored in metadata.

## Usage

``` r
filter_non_normal(exposomicset, p_thresh = 0.05)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure and omics data.

- p_thresh:

  A numeric value specifying the p-value threshold for normality.
  Variables with `p.value < p_thresh` are removed. Default is `0.05`.

## Value

A `MultiAssayExperiment` object with non-normal exposure variables
removed.

## Details

The function identifies exposure variables that fail a normality test
using `metadata(exposomicset)$transformation$norm_df`.

- Exposure variables with `p.value < p_thresh` are removed from
  `colData(exposomicset)`.

- The same filtering is applied to `colData` in each experiment within
  `experiments(exposomicset)`.

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

# Test for normality
mae <- mae |>
    run_normality_check() |>
    transform_exposure(exposure_cols = c("age", "bmi", "exposure_pm25"))
#> Checking Normality Using Shapiro-Wilk Test
#> 4 Exposure Variables are Normally Distributed
#> 0 Exposure Variables are NOT Normally Distributed
#> Applying the boxcox_best transformation.

# Remove non-normal variables
mae_filtered <- mae |>
    filter_non_normal()
#> Filtering out 0 non-normal exposure variables.
```
