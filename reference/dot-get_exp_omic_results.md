# Get Filtered Exposure-Omics Association Results

Internal helper that extracts the exposure-omics results and keeps rows
passing both the p-value and effect size thresholds.

## Usage

``` r
.get_exp_omic_results(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  estimate_thresh = 0
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing association results.

- pval_col:

  Character. Name of the column used for p-value filtering. Defaults to
  `"p_adjust"`.

- pval_thresh:

  Numeric. Significance threshold applied to `pval_col`. Rows with
  values below this threshold are retained. Defaults to `0.05`.

- estimate_thresh:

  Numeric. Minimum absolute effect size. Rows with `abs(estimate)`
  greater than or equal to this value are retained. Defaults to `0` (no
  filtering on estimate).

## Value

A tibble of filtered association results.
