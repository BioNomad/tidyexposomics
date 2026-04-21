# Plot Sample Clusters

Generates a heatmap of sample clustering results and summarizes sample
group assignments.

## Usage

``` r
plot_sample_clusters(exposomicset, exposure_cols = NULL)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing sample clustering results
  in `metadata(exposomicset)$sample_clustering`.

- exposure_cols:

  A character vector specifying columns from `colData` to include in the
  summary. Default is `NULL`, which includes all available columns.

## Value

A `ComplexHeatmap` plot displaying sample clustering results.

## Details

This function:

- Extracts sample cluster assignments from
  `metadata(exposomicset)$sample_clustering`.

- Merges cluster labels with `colData(exposomicset)`.

- Plots the heatmap stored in
  `metadata(exposomicset)$sample_clustering$heatmap`.

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

# determine sample clusters
mae <- run_cluster_samples(
    exposomicset = mae,
    exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
    clustering_approach = "diana"
)
#> Starting clustering analysis...
#> Optimal number of clusters for samples: 11

# plot sample clusters
sample_cluster_p <- mae |>
    plot_sample_clusters(
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )
#> tidyHeatmap says: If you use tidyHeatmap for scientific research, please cite: Mangiola, S. and Papenfuss, A.T., 2020. 'tidyHeatmap: an R package for modular heatmap production based on tidy principles.' Journal of Open Source Software. doi:10.21105/joss.02472.
#> This message is displayed once per session.
```
