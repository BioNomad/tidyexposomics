# Volcano Plot of Differential Abundance

Generates a **volcano plot** to visualize differential abundance results
across one or more omics layers.

## Usage

``` r
plot_volcano(
  exposomicset,
  pval_col = "adj.P.Val",
  pval_thresh = 0.05,
  logFC_col = "logFC",
  logFC_thresh = log2(1.5),
  plot_n_sig = TRUE,
  top_n_label = NULL,
  features_to_label = NULL,
  feature_col = "feature",
  xlab = expression(Log[2] * "FC"),
  ylab = expression(-Log[10] * "P"),
  title = "Volcano Plot of Differential Abundance",
  nrow = 2
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing differential abundance
  results in `metadata(exposomicset)$differential_abundance`.

- pval_col:

  A character string specifying the column containing p-values. Default
  is `"adj.P.Val"`.

- pval_thresh:

  A numeric threshold for significance. Features with p-values below
  this are considered significant. Default is `0.05`.

- logFC_col:

  A character string specifying the column for log fold changes. Default
  is `"logFC"`.

- logFC_thresh:

  A numeric threshold for absolute log fold change significance. Default
  is `log2(1.5)`.

- plot_n_sig:

  Logical; if `TRUE`, appends the number of significant features to
  facet titles. Default is `TRUE`.

- top_n_label:

  Optional integer. If provided, the top `n` most significant features
  per assay will be labeled on the plot.

- features_to_label:

  Optional character vector. Specific features to label regardless of
  significance.

- feature_col:

  A character string naming the feature ID column to use for labeling.
  Default is `"feature"`.

- xlab:

  Label for the x-axis. Default is `expression(Log[2]*"FC")`.

- ylab:

  Label for the y-axis. Default is `expression(-Log[10]*"P")`.

- title:

  Plot title. Default is `"Volcano Plot of Differential Abundance"`.

- nrow:

  Number of rows in the `facet_wrap()` layout. Default is `2`.

## Value

A `ggplot2` object representing the volcano plot.

## Details

The function:

- Extracts differential abundance results from
  `metadata(exposomicset)$differential_abundance`.

- Assigns each feature a direction of change: **Upregulated**,
  **Downregulated**, or **Not-Significant**.

- Uses `logFC_thresh` and `pval_thresh` to define thresholds.

- Adds dashed lines to indicate cutoffs for fold change and
  significance.

- Uses `facet_wrap()` to display each assay (`exp_name`) separately.

- Optionally labels the most significant features or user-defined ones.

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

# perform differential abundance analysis
mae <- run_differential_abundance(
    exposomicset = mae,
    formula = ~ smoker + sex,
    abundance_col = "counts",
    method = "limma_voom",
    action = "add"
)
#> Running differential abundance testing.
#> Processing assay: mRNA
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Processing assay: proteomics
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Differential abundance testing completed.

# create the volcano plot
volcano_p <- mae |>
    plot_volcano()
```
