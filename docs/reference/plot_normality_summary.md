# Plot Normality Summary of Exposure Variables

Generates a bar plot summarizing the number of exposure variables that
pass or fail normality tests (e.g., Shapiro-Wilk) before or after
transformation.

## Usage

``` r
plot_normality_summary(exposomicset, transformed = FALSE)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with quality control metadata.

- transformed:

  Logical; if `TRUE`, use results after transformation. Default is
  `FALSE`.

## Value

A `ggplot` object summarizing the number of exposures classified as
normal or not normal.

## Details

This function assumes that
[`run_normality_check()`](https://BioNomad.github.io/tidyexposomics/reference/run_normality_check.md)
has been executed and that the results are stored in
`metadata(exposomicset)$quality_control$normality`. If
`transformed = TRUE`, the function will instead plot the transformation
summary stored in
`metadata(exposomicset)$quality_control$transformation$norm_summary`,
which is populated by
[`transform_exposure()`](https://BioNomad.github.io/tidyexposomics/reference/transform_exposure.md).

The plot includes both bar heights and overlaid line segments to
reinforce the counts.

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

# plot the normality summary
norm_p <- mae |>
    plot_normality_summary()
```
