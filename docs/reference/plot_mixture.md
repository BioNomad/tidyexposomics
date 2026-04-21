# Plot Mixture Analysis Results

Create visualizations for mixture analysis results from qgcomp, WQS, or
BKMR.

## Usage

``` r
plot_mixture(
  exposomicset,
  method = c("qgcomp", "wqs", "bkmr"),
  plot_type = NULL,
  threshold = TRUE,
  top_n = NULL,
  ...
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with mixture analysis results, or a
  mixture analysis result object directly.

- method:

  Which method's results to plot: `"qgcomp"`, `"wqs"`, or `"bkmr"`.

- plot_type:

  Type of plot to generate. Options depend on method:

  For qgcomp: `"weights"` (default), `"effect"`

  For WQS: `"weights"` (default), `"effect"`

  For BKMR: `"pips"` (default), `"univariate"`, `"overall"`,
  `"interaction"`

- threshold:

  For weight plots, show reference line at 1/n_exposures. Default is
  `TRUE`.

- top_n:

  For weight plots, only show top n exposures by weight. Default is
  `NULL` (show all).

- ...:

  Additional arguments passed to underlying plot functions.

## Value

A ggplot2 object.

## Examples

``` r
mae <- make_example_data(n_samples = 100, return_mae = TRUE)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

mae <- mae |>
    run_mixture_analysis(
        outcome = "outcome_bmi",
        exposures = c("exposure_pm25", "exposure_no2"),
        method = "qgcomp",
        n_boot = 0
    )
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car 
#> Error in run_mixture_analysis(mae, outcome = "outcome_bmi", exposures = c("exposure_pm25",     "exposure_no2"), method = "qgcomp", n_boot = 0): Variables not found in colData: outcome_bmi

plot_mixture(mae, method = "qgcomp", plot_type = "weights")
#> Error in plot_mixture(mae, method = "qgcomp", plot_type = "weights"): No qgcomp results found. Run run_mixture_analysis() first.
```
