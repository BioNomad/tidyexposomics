# Transform Exposure Variables for Normality

Applies a transformation to selected numeric exposure variables in the
`colData` of a `MultiAssayExperiment` to improve their normality (e.g.,
log, Box-Cox, sqrt). Transformation results and normality statistics are
stored in metadata for tracking.

## Usage

``` r
transform_exposure(
  exposomicset,
  exposure_cols = NULL,
  transform_method = "boxcox_best"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure variables in
  `colData`.

- exposure_cols:

  Optional character vector of exposure variable names to transform. If
  `NULL`, uses exposures found in
  `metadata(exposomicset)$quality_control$normality$norm_df$exposure`.

- transform_method:

  Character. Transformation method to apply. Options:

  - `"none"`: no transformation

  - `"log2"`: log base 2 transformation

  - `"x_1_3"`: cube-root transformation

  - `"sqrt"`: square-root transformation

  - `"boxcox_best"`: data-driven Box-Cox approximation with heuristic
    labeling

## Value

A `MultiAssayExperiment` object with transformed exposures in `colData`,
and transformation details stored in:

- `metadata(exposomicset)$quality_control$transformation$norm_df`:
  Shapiro-Wilk test results

- `metadata(exposomicset)$quality_control$transformation$norm_summary`:
  Summary of normality

- `metadata(exposomicset)$codebook`: Updated with transformation info
  per variable

- `metadata(exposomicset)$summary$steps`: Updated with step record

## Details

For `transform_method = "boxcox_best"`, the function automatically
shifts values to be strictly positive and chooses from a discrete set of
transformations (e.g., `1/x`, `log(x)`, `x^2`) based on estimated
Box-Cox lambda. Each variable may receive a different transformation.

## See also

[`boxcox`](https://rdrr.io/pkg/MASS/man/boxcox.html),
[`shapiro.test`](https://rdrr.io/r/stats/shapiro.test.html)

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
    transform_exposure(
        exposure_cols = c("age", "bmi", "exposure_pm25"),
        transform_method = "boxcox_best"
    )
#> Checking Normality Using Shapiro-Wilk Test
#> 4 Exposure Variables are Normally Distributed
#> 0 Exposure Variables are NOT Normally Distributed
#> Applying the boxcox_best transformation.
```
