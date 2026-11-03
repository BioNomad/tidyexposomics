# Run functional enrichment on omics features

Runs over-representation analysis with fenr on features selected from an
`MultiAssayExperiment`-based exposomicset. Features are mapped to gene
identifiers, tested against a per-experiment universe (all features in
that experiment's `rowData`), p-values are adjusted across all tests,
and enriched terms are clustered by gene overlap (Jaccard similarity).

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
  action = c("add", "get")
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` (exposomicset) object with the upstream
  analysis results stored in its metadata.

- feature_type:

  Character. Which features to test:

  `"degs"`

  :   Differentially abundant features passing the p-value and log
      fold-change thresholds.

  `"degs_robust"`

  :   Features with a stability score above `score_threshold` from the
      sensitivity analysis.

  `"omics"`

  :   Features supplied via `variable_map`.

  `"factor_features"`

  :   Features loading on multi-omics integration factors (see
      `factor_type`).

  `"degs_cor"`, `"omics_cor"`, `"factor_features_cor"`

  :   Features from the corresponding correlation results, tested per
      `category`.

- score_col:

  Character. Column in the feature stability table used for
  `"degs_robust"`. Default `"stability_score"`.

- score_threshold:

  Numeric or `NULL`. Minimum stability score. If `NULL`, the threshold
  stored in the sensitivity analysis is used.

- variable_map:

  A data frame with a `variable` column (and `exp_name`, `category`
  where applicable) mapping variables to features. Required when
  `feature_type` is `"omics"`.

- factor_type:

  Character. Which factor feature set to use for `"factor_features"`:
  `"common_top_factor_features"` or `"top_factor_features"`.

- feature_col:

  Character. Column in the feature metadata (`rowData`) holding the
  identifier used for enrichment, e.g. gene symbols. Default
  `"feature"`.

- deg_pval_col:

  Character. Column with adjusted p-values in the differential abundance
  table. Default `"adj.P.Val"`.

- deg_pval_threshold:

  Numeric. P-value cutoff for calling DEGs (strictly less than). Default
  `0.05`.

- deg_logfc_col:

  Character. Column with log fold-changes in the differential abundance
  table. Default `"logFC"`.

- deg_logfc_threshold:

  Numeric. Absolute log fold-change cutoff for calling DEGs (strictly
  greater than). Default `log2(1.5)`.

- db:

  Character. Annotation database: `"GO"`, `"KEGG"`, or `"Reactome"`.

- species:

  Character or `NULL`. Species designation passed to fenr. Required for
  `db = "GO"`, and used for Reactome.

- fenr_col:

  Character. Name of the feature column in the fenr term mapping to
  match against. Default `"gene_symbol"`.

- padj_method:

  Character. Method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"fdr"`.

- pval_thresh:

  Numeric. Adjusted p-value cutoff for retained terms. Default `0.1`.

- min_set, max_set:

  Numeric. Terms are kept only if the number of selected features
  annotated to them is greater than `min_set` and less than `max_set`.
  Defaults `3` and `800`.

- clustering_approach:

  Character. Clustering method used to group enriched terms from the
  Jaccard similarity matrix. Default `"diana"`.

- action:

  Character. `"add"` stores the results in the object's metadata and
  returns the object; `"get"` returns the enrichment results directly.

## Value

If `action = "get"`, a tibble of enrichment results with `exp_name`,
adjusted p-values (`padj`), a `go_group` cluster assignment and, for
`_cor` feature types, `category`. If `action = "add"`, the input
`exposomicset` with results stored in
`metadata(exposomicset)$enrichment[[feature_type]]` and the run
parameters recorded in `metadata(exposomicset)$summary$steps`.

## Details

Term data are downloaded only when there are features to test. For
`db = "GO"`, fenr's default annotation URL is overridden and
[`fenr::fetch_go()`](https://rdrr.io/pkg/fenr/man/fetch_go.html) is
bypassed, because both rely on endpoints that are no longer available.
Fewer than two selected genes in an experiment/category skips that test;
fewer than two enriched terms skips clustering and sets `go_group` to
`NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Enrichment of DEGs against GO (biological process terms)
mae <- run_enrichment(
    exposomicset = mae,
    feature_type = "degs",
    db = "GO",
    species = "hsapiens"
)

# Return results without modifying the object
res <- run_enrichment(
    exposomicset = mae,
    feature_type = "degs_robust",
    db = "KEGG",
    action = "get"
)
} # }
```
