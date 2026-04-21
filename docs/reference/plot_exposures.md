# Plot Exposure Distributions by Category or Group

Visualizes exposure variable distributions using **boxplots** or **ridge
plots**.

## Usage

``` r
plot_exposures(
  exposomicset,
  exposure_cat = "all",
  exposure_cols = NULL,
  group_by = NULL,
  plot_type = "boxplot",
  alpha = 0.3,
  panel_sizes = rep(1, 100),
  title = "Exposure Levels by Category",
  xlab = "",
  ylab = "",
  facet_cols = NULL,
  group_cols = NULL,
  box_width = 0.1,
  fill_lab = ""
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure data.

- exposure_cat:

  A character string or vector specifying exposure category names (from
  `codebook$category`) to include. Use `"all"` to include all exposures.

- exposure_cols:

  Optional character vector specifying exact exposure variables to plot.

- group_by:

  A string specifying the column in `colData(exposomicset)` used to fill
  the plot (e.g., `"sex"`). Defaults to `NULL`, in which case exposures
  are colored by `category`.

- plot_type:

  Type of plot: `"boxplot"` (default) or `"ridge"`.

- alpha:

  Transparency level for background facet color strips. Default is
  `0.5`.

- panel_sizes:

  A numeric vector passed to
  [`ggh4x::force_panelsizes()`](https://teunbrand.github.io/ggh4x/reference/force_panelsizes.html)
  for controlling facet widths or heights.

- title:

  Plot title. Default is `"Exposure Levels by Category"`.

- xlab:

  X-axis label. Default is an empty string.

- ylab:

  Y-axis label. Default is an empty string.

- facet_cols:

  Optional vector of colors to use as background for facet categories.
  If `NULL`, a default palette is used.

- group_cols:

  Optional named vector of colors for `group_by` levels. If `NULL`, a
  default palette is used.

- box_width:

  A numeric value specifying the width of the boxplots. Only used when
  `plot_type = "boxplot"`. Default is `0.1`.

- fill_lab:

  Legend title for the fill aesthetic (e.g., `"Sex"` or
  `"Exposure Group"`). Default is `""`.

## Value

A `ggplot2` object showing exposure distributions, optionally grouped.

## Details

This function:

- Filters exposure data based on category or selected columns.

- Merges variable metadata from `metadata(exposomicset)$codebook`.

- Supports either **boxplot** (vertical distributions per variable) or
  **ridgeplot** (horizontal density plots per variable).

- If `group_by` is specified, that variable defines the plot fill color;
  otherwise, the fill is based on exposure `category`.

- Facets by `category` using
  [`ggh4x::facet_grid2()`](https://teunbrand.github.io/ggh4x/reference/facet_grid2.html)
  with color-coded strip backgrounds.

- The `box_width` argument controls the width of the boxplots when
  `plot_type = "boxplot"`.

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

# plot exposure data
exposure_plot <- mae |>
    plot_exposures(
        exposure_cols = c("exposure_pm25", "exposure_no2"),
        box_width = 0.2
    )
```
