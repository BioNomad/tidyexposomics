# Plot an Association Model Comparison

Visualise the output of
[`compare_association()`](https://BioNomad.github.io/tidyexposomics/reference/compare_association.md)
either as a dodged forest plot (one colour per model) or as a scatter of
each model's estimates against the reference model's.

## Usage

``` r
plot_association_comparison(
  comparison,
  type = c("forest", "scatter"),
  colors = "default",
  title = "Association Comparison"
)
```

## Arguments

- comparison:

  Output of
  [`compare_association()`](https://BioNomad.github.io/tidyexposomics/reference/compare_association.md).

- type:

  `"forest"` (default) or `"scatter"`. The scatter plots every
  non-reference model against the reference, with the identity line;
  points on the line are unchanged between models.

- colors:

  character vector of colors, used to color the points.

- title:

  Plot title.

## Value

A `ggplot2` object.

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

mae |>
    compare_association(source = "exposures") |>
    plot_association_comparison()
```
