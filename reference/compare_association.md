# Compare Association Models

Line up the results of two or more stored models (e.g. crude vs
adjusted, or female vs male strata) feature by feature.

## Usage

``` r
compare_association(
  exposomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  models = NULL,
  reference = NULL,
  effect = "main",
  terms = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` with results from
  [`run_association()`](https://BioNomad.github.io/tidyexposomics/reference/run_association.md).

- source:

  Association source: `"omics"`, `"exposures"`, `"factors"`, or
  `"go_pcs"`.

- models:

  Character vector of stored model names to compare. Default `NULL`
  compares every model stored for the source.

- reference:

  Model to use as the reference for the derived columns. Defaults to the
  first of `models`.

- effect:

  Which effect to compare: `"main"` (default) or an interacting
  coefficient name such as `"sexMale"`.

- terms:

  Optional character vector of features to keep.

## Value

A tibble with one row per feature and model: `term`, `model`,
`relation`, `estimate`, `std.error`, `p.value`, any of `p_adjust`, `n`,
`n_events`, then `ref_estimate`, `ref_se`, `reference`, `pct_change`,
`z_diff`, and `p_diff`.

## Details

The `relation` column says how each model relates to the reference,
which determines which derived statistics are filled in:

- `"same_sample"`: same samples, different model (e.g. crude vs
  adjusted). `pct_change` is the change-in-estimate,
  `100 * (estimate - reference estimate) / reference estimate`; negative
  values mean the effect shrank in magnitude.

- `"sibling_strata"`: different strata from the same `strata` run.
  `z_diff` and `p_diff` test the difference between the two estimates,
  `(b1 - b2) / sqrt(se1^2 + se2^2)`, which assumes the strata are
  independent. Comparisons are pairwise against the reference.

- `"unrelated"`: no derived statistics.

## Examples

``` r
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
) |>
    run_association(
        formula = list(
            crude = smoker ~ 1,
            adjusted = smoker ~ age
        ),
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial"
    )
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.
#> Running GLMs.

compare_association(
    mae,
    source = "exposures"
)
#> # A tibble: 4 × 15
#>   term         model relation estimate std.error p.value p_adjust     n n_events
#>   <chr>        <fct> <chr>       <dbl>     <dbl>   <dbl>    <dbl> <int>    <dbl>
#> 1 exposure_no2 crude referen…   -0.577     0.698   0.409    0.409    10        5
#> 2 exposure_no2 adju… same_sa…   -0.650     0.738   0.378    0.393    10        5
#> 3 exposure_pm… crude referen…   -0.629     0.717   0.380    0.409    10        5
#> 4 exposure_pm… adju… same_sa…   -0.619     0.724   0.393    0.393    10        5
#> # ℹ 6 more variables: ref_estimate <dbl>, ref_se <dbl>, reference <chr>,
#> #   pct_change <dbl>, z_diff <dbl>, p_diff <dbl>
```
