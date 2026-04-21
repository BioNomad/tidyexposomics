# Run Mixture Analysis

Perform mixture analysis to assess joint effects of multiple exposures
on an outcome using quantile g-computation (qgcomp), weighted quantile
sum (WQS), or Bayesian kernel machine regression (BKMR).

## Usage

``` r
run_mixture_analysis(
  exposomicset,
  outcome,
  exposures = NULL,
  covariates = NULL,
  method = c("qgcomp", "wqs", "bkmr"),
  family = "gaussian",
  n_quantiles = 4,
  n_boot = 200,
  repeat_holdout = 1,
  validation = 0.6,
  direction = "both",
  lambda = 0,
  n_iter = 10000,
  varsel = TRUE,
  degree = 1,
  id = NULL,
  weights = NULL,
  seed = 123,
  parallel = FALSE,
  keep_fit = FALSE,
  action = "add",
  ...
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object.

- outcome:

  The outcome variable name (must be in `colData`).

- exposures:

  Character vector of exposure variable names. If `NULL`, uses all
  numeric variables from the codebook.

- covariates:

  Optional character vector of covariate names.

- method:

  One of `"qgcomp"`, `"wqs"`, or `"bkmr"`.

- family:

  `"gaussian"`, `"binomial"`, `"poisson"`, or `"cox"`.

- n_quantiles:

  Number of quantiles. Default is 4.

- n_boot:

  Number of bootstrap samples for qgcomp/WQS. Default is 200. For
  qgcomp, set to 0 for no bootstrap (faster, uses Delta method for CI).

- repeat_holdout:

  Number of repeated holdout validations for WQS. Default is 1.

- validation:

  Proportion of data for validation in WQS. Default is 0.6.

- direction:

  For WQS: `"positive"`, `"negative"`, or `"both"`.

- lambda:

  Penalization parameter for WQS. Default is 0 (no penalty).

- n_iter:

  Number of iterations for BKMR. Default is 10000.

- varsel:

  Logical. Whether to perform variable selection in BKMR. Default is
  `TRUE`.

- degree:

  Polynomial degree for qgcomp MSM. Default is 1 (linear).

- id:

  Optional cluster ID variable name for robust standard errors.

- weights:

  Optional sampling weights variable name.

- seed:

  Random seed for reproducibility.

- parallel:

  Logical. Whether to use parallel processing for bootstrap. Default is
  `FALSE`.

- keep_fit:

  Logical. Whether to keep the full model fit object. Default is
  `FALSE`. Set to `TRUE` if you need access to the underlying model for
  diagnostics or plotting.

- action:

  If `"add"` (default), saves results to metadata. if `"get"`, returns
  results as list.

- ...:

  Additional arguments passed to the underlying method.

## Value

If `action = "add"`, returns updated `MultiAssayExperiment`. Otherwise,
returns a list containing:

For all methods:

- `weights`: tibble of exposure weights

- `mixture_effect`: tibble with overall mixture effect estimate

- `method`, `outcome`, `exposures`, `covariates`: analysis metadata

Method-specific outputs:

- qgcomp: `partial_effects` tibble

- wqs: `model_summary` tibble

- bkmr: `pips` tibble (posterior inclusion probabilities)

## Details

Methods available:

- `qgcomp`: Quantile g-computation via `qgcomp` package. This is a fast
  method and allows negative weights, estimates effects in both
  directions simultaneously. Recommended as the default method. Supports
  clustering via `id` parameter and sampling weights via `weights`
  parameter.

- `wqs`: Weighted Quantile Sum regression via `gWQS` package.
  Constrained weights (need to be all positive or all negative), this is
  good for directional hypotheses.

- `bkmr`: Bayesian Kernel Machine Regression via `bkmr` package. This is
  a flexible nonparametric approach, which can detect nonlinearity and
  interactions, but is a slower method.

## Examples

``` r
mae <- make_example_data(n_samples = 100, return_mae = TRUE)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# Quantile g-computation (fast, recommended)
mae |>
    run_mixture_analysis(
        outcome = "outcome_bmi",
        exposures = c("exposure_pm25", "exposure_no2"),
        method = "qgcomp"
    )
#> Error in run_mixture_analysis(mae, outcome = "outcome_bmi", exposures = c("exposure_pm25",     "exposure_no2"), method = "qgcomp"): Variables not found in colData: outcome_bmi

# qgcomp without bootstrap (even faster)
mae |>
    run_mixture_analysis(
        outcome = "outcome_bmi",
        exposures = c("exposure_pm25", "exposure_no2"),
        method = "qgcomp",
        n_boot = 0
    )
#> Error in run_mixture_analysis(mae, outcome = "outcome_bmi", exposures = c("exposure_pm25",     "exposure_no2"), method = "qgcomp", n_boot = 0): Variables not found in colData: outcome_bmi

# WQS with both directions
mae |>
    run_mixture_analysis(
        outcome = "outcome_bmi",
        exposures = c("exposure_pm25", "exposure_no2"),
        method = "wqs",
        direction = "both"
    )
#> Error in run_mixture_analysis(mae, outcome = "outcome_bmi", exposures = c("exposure_pm25",     "exposure_no2"), method = "wqs", direction = "both"): Variables not found in colData: outcome_bmi
```
