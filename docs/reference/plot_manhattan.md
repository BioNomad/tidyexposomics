# Plot a Manhattan-style ExWAS summary across omics categories

This function generates a multi-faceted Manhattan plot from the results
of `associate_all_outcome()`, visualizing the significance of
associations across omics features, grouped by category. Significant
features can be highlighted and labeled, and strip backgrounds can be
colored per facet.

## Usage

``` r
plot_manhattan(
  exposomicset,
  pval_thresh = 0.05,
  feature_col = "term",
  alpha = 0.5,
  min_per_cat = 1,
  vars_to_label = NULL,
  sig_color = "magenta2",
  non_sig_cols = c("grey25", "grey75"),
  pval_thresh_line_col = "grey25",
  panel_sizes = c(1, 1, 1, 1, 1),
  linetype = "dashed",
  facet_cols = NULL,
  label_size = 3.5,
  facet_angle = 90,
  facet_text_face = "bold.italic"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object that has already been processed by
  `associate_all_outcome()`.

- pval_thresh:

  Numeric threshold for significance (default = 0.05).

- feature_col:

  A character string indicating the column name to use for feature
  labeling and highlighting (e.g., `"term"` or `"feature"`). Default is
  `"term"`.

- alpha:

  Transparency applied to facet strip colors (default = 0.5).

- min_per_cat:

  Minimum number of features per category to be shown (default = 1).

- vars_to_label:

  Optional character vector of variable names to label explicitly,
  matched against the `feature_col` column.

- sig_color:

  Color used for significant points (default = `"magenta2"`).

- non_sig_cols:

  Character vector of alternating colors for non-significant points
  (default = `c("grey25", "grey75")`).

- pval_thresh_line_col:

  Color of the horizontal significance threshold line (default =
  `"grey25"`).

- panel_sizes:

  Numeric vector passed to
  [`ggh4x::force_panelsizes()`](https://teunbrand.github.io/ggh4x/reference/force_panelsizes.html)
  to control panel widths (default = `c(1,1,1,1,1)`).

- linetype:

  Line type for the horizontal threshold (default = `"dashed"`).

- facet_cols:

  Optional vector of colors to use for facet strip backgrounds.

- label_size:

  Numeric size of the feature label text (default = 3.5).

- facet_angle:

  Angle (in degrees) for strip text rotation (default = 90).

- facet_text_face:

  Font face for facet strip labels (default = `"bold.italic"`).

## Value

A `ggplot` object showing the Manhattan-style faceted plot.

## Details

- This function expects `associate_all_outcome()` to have been run
  first.

- Facets represent omics categories, and points represent features.

- Points below the significance threshold are colored using
  `non_sig_cols`, while significant ones are colored with `sig_color`
  and optionally labeled.

- Uses `ggrepel` to avoid overlapping labels and `ggh4x` for enhanced
  faceting.

- The `feature_col` argument allows customization of which column is
  used to label or identify features, enabling compatibility with
  different result formats.

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
        source = "omics",
        top_n = 20,
        feature_set = c("exposure_pm25", "exposure_no2"),
        outcome = "smoker",
        covariates = c("age"),
        family = "binomial"
    )
#> Log2-Transforming each assay in MultiAssayExperiment.
#> Scaling each assay in MultiAssayExperiment.
#> Running GLMs.
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

# create the manhattan plot
manhattan_p <- mae |>
    plot_manhattan(feature_col = "term")
```
