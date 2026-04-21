# Plot Significant Omics Associations per Exposure

Internal helper that counts and plots the number of significant omics
associations for each individual exposure.

## Usage

``` r
.plot_exp_omic_exposure(
  exposomicset,
  pval_col = "p_adjust",
  pval_thresh = 0.05
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

## Value

A `ggplot` object.
