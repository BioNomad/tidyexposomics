# Summarize and Visualize Analysis Pipeline Steps

This function prints and visualizes the analysis steps stored in the
metadata of a `MultiAssayExperiment` object. The steps are optionally
printed to the console as a numbered list and can be rendered as a
left-to-right Mermaid flowchart. The flowchart connects steps with
arrows and includes step notes if requested.

## Usage

``` r
run_pipeline_summary(
  exposomicset,
  show_index = TRUE,
  console_print = TRUE,
  diagram_print = FALSE,
  include_notes = TRUE
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object that contains a "summary" entry in its
  metadata, which includes a list named `steps`.

- show_index:

  Logical, default `TRUE`. If `TRUE`, prefixes each step with its index.

- console_print:

  Logical, default `TRUE`. If `TRUE`, prints the step list to the
  console.

- diagram_print:

  Logical, default `FALSE`. If `TRUE`, renders a Mermaid diagram of the
  steps.

- include_notes:

  Logical, default `TRUE`. If `TRUE`, appends any "notes" associated
  with each step to the label.

## Value

No return value. This function is called for its side effects: console
output and/or diagram rendering.

## Details

The Mermaid flowchart is rendered left-to-right and connects each step
in sequence. Each node is labeled using the step name and, optionally,
any attached notes.

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

# Test for normality
mae <- mae |>
    run_normality_check() |>
    transform_exposure(exposure_cols = c("age", "bmi", "exposure_pm25"))
#> Checking Normality Using Shapiro-Wilk Test
#> 3 Exposure Variables are Normally Distributed
#> 1 Exposure Variables are NOT Normally Distributed
#> Applying the boxcox_best transformation.

# Run the pipeline summary
run_pipeline_summary(mae)
#> 1. run_normality_check - Assessed normality of 4 numeric exposure variables. 3 were normally distributed (p > 0.05), 1 were not.
#> 2. transform_exposure - Applied 'boxcox_best' transformation to 3 exposure variables. 2 passed normality (Shapiro-Wilk p > 0.05, 66.7%).
```
