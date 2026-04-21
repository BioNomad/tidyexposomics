# Plot Missing Data Across Exposure and Omic Layers

Visualizes missing data patterns in a `MultiAssayExperiment` object
using summary bar plots or feature-level lollipop plots.

## Usage

``` r
plot_missing(
  exposomicset,
  threshold = 5,
  plot_type = c("summary", "lollipop"),
  layers = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure and omics assays.
  Missing data is inferred directly from the assays.

- threshold:

  Numeric. The percentage threshold (0-100) above which features are
  counted as missing in the summary plot. Default is `5`.

- plot_type:

  Character. Type of plot to generate. Either `"summary"` for a bar plot
  showing number of features above the missing threshold, or
  `"lollipop"` for a per-feature lollipop plot with layer annotations.
  Default is `"summary"`.

- layers:

  Optional character vector. If specified, filters the plot to include
  only selected layers (e.g., `"Exposure"`, `"Transcriptome"`).

## Value

A `ggplot` or `patchwork` object depending on the selected `plot_type`.

## Details

The function calculates missing data per feature (or variable) across
all assays (including exposure variables) and generates:

- **Summary plot (`plot_type = "summary`)**: A bar plot showing the
  number of variables in each assay exceeding the specified missingness
  threshold.

- **Lollipop plot (`plot_type = "lollipop`)**: A feature-level plot
  where each feature's percent missingness is shown, along with a
  color-coded tile on the side indicating its layer of origin.

The tile colors in the lollipop plot match the experiment colors used in
other visualizations (e.g., via `scale_color_tidy_exp()`).

## Examples

``` r
#' # Create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# Introduce some missingness
MultiAssayExperiment::colData(mae)$exposure_pm25[sample(1:20, 5)] <- NA

# Summary bar plot of missing data
summary_p <- plot_missing(
    mae,
    threshold = 10,
    plot_type = "summary"
)

# Lollipop plot for all features with any missingness
lollipop_p <- plot_missing(
    mae,
    plot_type = "lollipop"
)
```
