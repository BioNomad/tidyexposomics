# Plot Association Results (Unified Forest Plot)

Generates a forest plot for association results from any source stored
in the metadata of a `MultiAssayExperiment` object. Supports faceting
and visual augmentation with R^2 tiles when available.

## Usage

``` r
plot_association(
  exposomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  terms = NULL,
  filter_col = "p.value",
  filter_thresh = 0.05,
  direction_filter = "all",
  add_r2_tile = TRUE,
  r2_col = "adj_r2",
  facet = FALSE,
  nrow = 1,
  subtitle = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing association results in
  metadata.

- source:

  Character string indicating the association source. One of `"omics"`,
  `"exposures"`, `"factors"`, or `"go_pcs"`.

- terms:

  Optional character vector of term names to subset the plot to. Default
  is `NULL` (include all).

- filter_col:

  Column used to assess statistical significance (default: `"p.value"`).

- filter_thresh:

  Numeric threshold applied to `filter_col` (default: `0.05`).

- direction_filter:

  Direction of associations to retain. One of `"all"` (default), `"up"`,
  or `"down"`.

- add_r2_tile:

  Logical; if `TRUE`, includes a tile plot for `r2_col` (default:
  `TRUE`).

- r2_col:

  Column used for coloring the tile plot (default: `"adj_r2"`).

- facet:

  Logical; if `TRUE` and `source == "go_pcs"`, apply nested faceting by
  experiment and GO cluster (default: `FALSE`).

- nrow:

  Integer; number of rows for facet layout if enabled (default: `1`).

- subtitle:

  Optional subtitle for the plot. If `NULL`, automatically generated
  from covariates used in the model.

## Value

A `ggplot2` object: either a single forest plot or a composite plot with
an R^2 tile strip.

## Details

This function visualizes effect size estimates and confidence intervals
from association analyses. It allows filtering by direction (`"up"` for
positive, `"down"` for negative) and by significance. For
`source = "go_pcs"`, it supports special formatting by splitting term
labels into nested facets.

The R^2 tile (if enabled) adds a side heatmap indicating model fit for
each association. This can be useful for model diagnostics.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# run association tests
mae <- mae |>
    run_association(
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        outcome = "smoker",
        covariates = c("age"),
        family = "binomial"
    )
#> Running GLMs.

assoc_plot <- mae |>
    plot_association(
        source = "exposures"
    )
```
