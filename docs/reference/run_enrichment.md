# Perform enrichment analysis on selected features from a exposomicset object

This function performs enrichment analysis using selected features
derived from differential expression, correlation analysis, or
multi-omics factor features across experiments in an `exposomicset`. It
supports multiple enrichment databases (e.g., GO, KEGG, Reactome),
applies FDR correction, and optionally clusters GO terms by Jaccard
overlap.

## Usage

``` r
run_enrichment(
  exposomicset,
  feature_type = c("degs", "degs_robust", "omics", "factor_features", "degs_cor",
    "omics_cor", "factor_features_cor"),
  score_col = "stability_score",
  score_threshold = NULL,
  variable_map = NULL,
  factor_type = c("common_top_factor_features", "top_factor_features"),
  feature_col = "feature",
  deg_pval_col = "adj.P.Val",
  deg_pval_threshold = 0.05,
  deg_logfc_col = "logFC",
  deg_logfc_threshold = log2(1.5),
  db = c("GO", "KEGG", "Reactome"),
  species = NULL,
  fenr_col = "gene_symbol",
  padj_method = "fdr",
  pval_thresh = 0.1,
  min_set = 3,
  max_set = 800,
  clustering_approach = "diana",
  action = "add"
)
```

## Arguments

- exposomicset:

  An `exposomicset` (a `MultiAssayExperiment` object with metadata)
  containing omics and metadata.

- feature_type:

  Character string indicating the feature source. One of `"degs"`,
  `"degs_robust"`, `"omics"`, `"factor_features"`, `"degs_cor"`,
  `"omics_cor"`, or `"factor_features_cor"`.

- score_col:

  Column name used for stability score filtering (only for
  `degs_robust`).

- score_threshold:

  Optional numeric threshold for filtering stability scores. If `NULL`,
  the default threshold stored in the metadata will be used.

- variable_map:

  A data frame with `exp_name` and `variable` columns, used when
  `feature_type = "omics"`.

- factor_type:

  Character string for selecting factor features:
  `"common_top_factor_features"` or `"top_factor_features"`.

- feature_col:

  The name of the feature column used to extract gene identifiers.

- deg_pval_col:

  Column name for adjusted p-values from DEG analysis.

- deg_pval_threshold:

  Threshold to select significant DEGs (default: 0.05).

- deg_logfc_col:

  Column name for log-fold changes from DEG analysis.

- deg_logfc_threshold:

  Threshold to select DEGs by absolute logFC (default: `log2(1.5)`).

- db:

  Enrichment database to use. One of `"GO"`, `"KEGG"`, `"Reactome"`.

- species:

  Species name (required for GO enrichment, e.g., `"Homo sapiens"`).
  Ignored for other databases.

- fenr_col:

  Column name for gene IDs used by `fenr` (e.g., `"gene_symbol"`).

- padj_method:

  Method for p-value adjustment (default: `"fdr"`).

- pval_thresh:

  Adjusted p-value threshold for filtering enriched terms (default:
  0.1).

- min_set:

  Minimum number of selected genes overlapping an enriched term
  (default: 3).

- max_set:

  Maximum number of selected genes overlapping an enriched term
  (default: 800).

- clustering_approach:

  Clustering method for GO term grouping. Defaults to `"diana"`.

- action:

  Either `"add"` to store results in the object's metadata or `"get"` to
  return results as a data frame.

## Value

If `action = "add"`, returns the modified `exposomicset` with enrichment
results added to metadata. If `action = "get"`, returns a `data.frame`
of enrichment results with GO term clusters (if applicable).

## Details

The function identifies selected features based on the chosen
`feature_type`, determines the gene universe for each experiment, and
performs enrichment analysis using the `fenr` package. Results are
adjusted for multiple testing and optionally clustered by gene set
overlap (for GO terms).

If `feature_type` includes correlation-based results (ending in `_cor`),
enrichment is performed for each exposure category separately.

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

# perform enrichment analysis
mae <- run_enrichment(
    exposomicset = mae,
    feature_type = "degs",
    feature_col = "symbol",
    species = "goa_human",
    deg_logfc_threshold = log2(1),
    deg_pval_col = "P.Value",
    deg_pval_threshold = 0.2,
    action = "add"
)
#> 
#> 
#> 
#> 
#> Determining Number of GO Term Clusters...
#> Optimal number of clusters for samples: 2
```
