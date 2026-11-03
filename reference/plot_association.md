# Plot Association Results

Generates a forest plot for association results from any source stored
in the metadata of a `MultiAssayExperiment` object. Supports faceting
and visual augmentation with R^2 tiles when available.

## Usage

``` r
plot_association(
  exposomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  model = NULL,
  terms = NULL,
  effect = "main",
  filter_col = "p.value",
  filter_thresh = 0.05,
  direction_filter = "all",
  add_r2_tile = TRUE,
  r2_col = "adj_r2",
  up_col = "#8E0152",
  down_col = "#006666",
  ns_col = "grey55",
  r2_low_col = "#e8e1eb",
  r2_high_col = "#3f1b4f",
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

- model:

  Name of the stored model to plot (see
  [`association_models()`](https://BioNomad.github.io/tidyexposomics/reference/association_models.md)).
  Can be left as `NULL` when only one model is stored for the source.

- terms:

  Optional character vector of term names to subset the plot to. Default
  is `NULL` (include all).

- effect:

  Character; which effect to plot. `"main"` (default) is the feature's
  main effect. For models with interactions, supply the name of the
  interacting coefficient (e.g. `"sexMale"`) to plot the interaction
  instead.

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

- up_col:

  Character; color used for significant positive associations (default:
  `"#8E0152"`).

- down_col:

  Character; color used for significant negative associations (default:
  `"#006666"`).

- ns_col:

  Character; color used for non-significant associations (default:
  `"grey55"`).

- r2_low_col:

  Character; color for the low end of the R^2 tile gradient (default:
  `"#e8e1eb"`). Only used when `add_r2_tile = TRUE`.

- r2_high_col:

  Character; color for the high end of the R^2 tile gradient (default:
  `"#3f1b4f"`). Only used when `add_r2_tile = TRUE`.

- facet:

  Logical; if `TRUE` and `source == "go_pcs"`, apply nested faceting by
  experiment and GO cluster (default: `FALSE`).

- nrow:

  Integer; number of rows for facet layout if enabled (default: `1`).

- subtitle:

  Optional subtitle for the plot. If `NULL`, automatically generated
  from the model formula and sample subset.

## Value

A `ggplot2` object: either a single forest plot or a composite plot with
an R^2 tile strip.

## Details

This function visualizes effect size estimates and confidence intervals
from association analyses. It allows filtering by direction (`"up"` for
positive, `"down"` for negative) and by significance.

The R^2 tile (if enabled) adds a side heatmap indicating model fit for
each association. This can be useful for model diagnostics.

To see several models together, use
[`compare_association()`](https://BioNomad.github.io/tidyexposomics/reference/compare_association.md)
and
[`plot_association_comparison()`](https://BioNomad.github.io/tidyexposomics/reference/plot_association_comparison.md).

## Examples

``` r
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

mae <- mae |>
    run_association(
        formula = smoker ~ age,
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial"
    )
#> Running GLMs.

assoc_plot <- mae |>
    plot_association(
        source = "exposures"
    )
```
