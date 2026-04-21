# Assess Normality of Exposure Variables

Performs Shapiro-Wilk tests to check the normality of numeric exposure
variables in `colData` of a `MultiAssayExperiment` object.

## Usage

``` r
run_normality_check(exposomicset, action = "add")
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure data in `colData`.

- action:

  A character string specifying whether to store (`"add"`) or return
  (`"get"`) the results. Default is `"add"`.

## Value

A `MultiAssayExperiment` object with normality results added to metadata
(if `action = "add"`) or a list with:

- norm_df:

  A data frame of Shapiro-Wilk test results for each exposure variable.

- norm_plot:

  A ggplot object showing the distribution of normal and non-normal
  exposures.

## Details

This function:

- Extracts **numeric, non-constant** exposure variables from `colData`.

- Runs **Shapiro-Wilk tests** to assess normality.

- Summarizes the number of normally and non-normally distributed
  exposures.

- Generates a bar plot visualizing the normality results.

- **Output Handling**:

  - `"add"`: Stores results in `metadata(exposomicset)$normality`.

  - `"get"`: Returns a list containing the normality test results and
    plot.

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
```
