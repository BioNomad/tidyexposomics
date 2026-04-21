# Plot PCA Results for Features and Samples

Generates PCA plots for both feature space and sample space, including
scatter plots and scree plots.

## Usage

``` r
plot_pca(
  exposomicset,
  feature_col = "#00a9b2",
  sample_col = "#8a4f77",
  sample_outlier_col = "firebrick"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing PCA results in
  `metadata(exposomicset)$pca`.

- feature_col:

  A character string specifying the color for the feature scree plot.
  Default is `"#00a9b2"`.

- sample_col:

  A character string specifying the color for the sample scree plot.
  Default is `"#8a4f77"`.

- sample_outlier_col:

  A character string specifying the color for sample outlier labels.
  Default is `"firebrick"`.

## Value

A combined `ggplot` object containing the four PCA plots.

## Details

This function creates four PCA visualizations:

- **Feature Space PCA Plot**: Colored by category (e.g., omics,
  exposure).

- **Feature Scree Plot**: Displays the variance explained by each
  principal component.

- **Sample Space PCA Plot**: Highlights outlier samples.

- **Sample Scree Plot**: Displays variance explained in the sample PCA.

Outliers are labeled based on `metadata(exposomicset)$pca$outliers`.

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

# run pca
mae <- mae |>
    run_pca()
#> Identifying common samples.
#> Subsetting exposure data.
#> Subsetting omics data.
#> Performing PCA on Feature Space.
#> Performing PCA on Sample Space.
#> No outliers detected.

# create the pca plot
pca_p <- mae |>
    plot_pca()
```
