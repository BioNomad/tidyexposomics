# Export tidyexposomics Results to Excel

Exports selected results stored in a `MultiAssayExperiment` object
created by the `tidyexposomics` pipeline to an Excel workbook. Users can
select which result types to include, and optionally add placeholder
sheets for missing data.

## Usage

``` r
extract_results_excel(
  exposomicset,
  file = "tidyexposomics_results.xlsx",
  result_types = c("correlation", "association", "mixture_analysis",
    "differential_analysis", "multiomics_integration", "network", "enrichment",
    "exposure_summary", "pipeline"),
  include_empty_tabs = FALSE
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with results stored in `@metadata`,
  typically created by the `tidyexposomics` pipeline.

- file:

  Character. Path to the output Excel file.

- result_types:

  Character vector specifying which result categories to export. Options
  include:

  - `"correlation"`: Correlation results.

  - `"association"`: Association results.

  - `"mixture_analysis"`: Mixture Analysis results.

  - `"differential_analysis"`: Differential abundance results, including
    sensitivity analysis if available.

  - `"multiomics_integration"`: Common top features contributing to
    latent factors from multi-omics integration.

  - `"network"`: Exposure impact metrics from network analyses.

  - `"enrichment"`: Enrichment results by omic and exposure category.

  - `"exposure_summary"`: Summary statistics for exposure variables.

  - `"pipeline"`: Overview of steps completed in the pipeline.

  Use `"all"` to export all of the above categories.

- include_empty_tabs:

  Logical. If `TRUE`, adds placeholder sheets for any missing result
  types. Default is `FALSE`.

## Value

An Excel file is written to the specified path. A message is printed
with the file location.

## Examples

``` r
# Create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# run correlation analysis
mae <- mae |>
    run_correlation(
        feature_type = "exposures",
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )

# file path of the output file
tmp <- tempfile(fileext = ".xlsx")

# extract the correlation results
extract_results_excel(
    exposomicset = mae,
    result_types = "correlation",
    file = tmp
)
#> Writing Correlation Results.
#> Results written to: C:\Users\Jason\AppData\Local\Temp\RtmpCm3inJ\file6b50286a2829.xlsx
```
