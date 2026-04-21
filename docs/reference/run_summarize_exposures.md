# Summarize Exposure Variables

Computes summary statistics for numeric exposure variables and
optionally stores the results in the `MultiAssayExperiment` metadata.

## Usage

``` r
run_summarize_exposures(exposomicset, exposure_cols = NULL, action = "add")
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure data in the sample
  metadata.

- exposure_cols:

  A character vector of exposure variable names to summarize. If `NULL`,
  all numeric exposure variables are included.

- action:

  A string specifying the action to take. Use `"add"` to attach the
  summary table to `metadata(exposomicset)` or `"get"` to return the
  summary table directly. Default is `"add"`.

## Value

A modified `MultiAssayExperiment` object (if `action = "add"`), or a
data frame of summary statistics (if `action = "get"`).

## Details

This function:

- Extracts sample-level exposure data using
  [`pivot_sample()`](https://BioNomad.github.io/tidyexposomics/reference/pivot_sample.md).

- Filters to user-specified exposures (`exposure_cols`) if provided.

- Computes descriptive statistics for each numeric variable:

  - Number of values (`n_values`)

  - Number of NAs (`n_na`)

  - Minimum, maximum, and range

  - Sum, median, mean

  - Standard error of the mean

  - 95% confidence interval of the mean

  - Variance, standard deviation

  - Coefficient of variation (`sd / mean`)

- Merges the result with variable metadata stored in
  `metadata(exposomicset)$codebook`.

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

# Summarize exposure data
exp_sum <- mae |>
    run_summarize_exposures(
        exposure_cols = c("age", "bmi", "exposure_pm25"),
        action = "get"
    )
```
