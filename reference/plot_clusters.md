# Plot Sample or Exposure Clusters

Generates a heatmap of either sample clustering or exposure clustering
results, split by cluster group.

## Usage

``` r
plot_clusters(
  exposomicset,
  cluster_type = c("sample", "exposure"),
  exposure_cols = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing clustering results in
  `metadata(exposomicset)$quality_control`.

- cluster_type:

  A character string specifying which clustering to plot: `"sample"`
  (from
  [`run_cluster_samples()`](https://BioNomad.github.io/tidyexposomics/reference/run_cluster_samples.md))
  or `"exposure"` (from
  [`run_cluster_exposures()`](https://BioNomad.github.io/tidyexposomics/reference/run_cluster_exposures.md)).
  Default is `"sample"`.

- exposure_cols:

  A character vector specifying columns from `colData` to include in the
  heatmap. Default is `NULL`, which includes all `colData` columns for
  `"sample"` and all clustered exposures for `"exposure"`.

## Value

A `ComplexHeatmap` plot displaying clustering results.

## Details

This function:

- Extracts cluster assignments from
  `metadata(exposomicset)$quality_control`.

- Merges cluster labels with `colData(exposomicset)` and the codebook.

- Plots a heatmap split by sample groups (columns) or exposure groups
  (rows), with annotation tiles for exposure category and cluster group.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 30,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

exposure_cols <- c("exposure_pm25", "exposure_no2", "age", "bmi")

# determine sample and exposure clusters
mae <- mae |>
    run_cluster_samples(
        exposure_cols = exposure_cols,
        clustering_approach = "diana"
    ) |>
    run_cluster_exposures(
        exposure_cols = exposure_cols,
        clustering_approach = "diana"
    )
#> Starting exposure clustering analysis...
#> Starting clustering analysis...
#> Optimal number of clusters for samples: 22
#> Optimal number of clusters for exposures: 1

# plot sample clusters
sample_cluster_p <- mae |>
    plot_clusters(
        cluster_type = "sample",
        exposure_cols = exposure_cols
    )
#> tidyHeatmap says: If you use tidyHeatmap for scientific research, please cite: Mangiola, S. and Papenfuss, A.T., 2020. 'tidyHeatmap: an R package for modular heatmap production based on tidy principles.' Journal of Open Source Software. doi:10.21105/joss.02472.
#> This message is displayed once per session.
#> Warning: `when()` was deprecated in purrr 1.0.0.
#> ℹ Please use `if` instead.
#> ℹ The deprecated feature was likely used in the tidyHeatmap package.
#>   Please report the issue at
#>   <https://github.com/stemangiola/tidyHeatmap/issues>.

# plot exposure clusters
exposure_cluster_p <- mae |>
    plot_clusters(
        cluster_type = "exposure",
        exposure_cols = exposure_cols
    )
#> Error in h(simpleError(msg, call)): error in evaluating the argument '.data' in selecting a method for function 'annotation_tile': error in evaluating the argument '.data' in selecting a method for function 'annotation_tile': tidyHeatmap says: the grouping exposure_group is not specific to row or columns. Maybe you just have one grouping.
```
