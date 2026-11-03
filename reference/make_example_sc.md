# Make an example single-cell dataset with multiple exposures

Creates a sample `SingleCellExperiment` object with sample-level
exposures (two continuous, one binary, plus age), cell types, and
negative binomial counts. A subset of genes truly responds to each
exposure, and some exposures act only in specific cell types. The true
effect sizes are stored in `rowData()` so results can be benchmarked.

## Usage

``` r
make_example_sc(
  n_samples = 12,
  cells_per_sample = 200,
  n_genes = 1000,
  n_celltypes = 4,
  prop_responsive = 0.08,
  effect_sd = 0.6,
  prop_celltype_responsive = 0.6,
  size = 2
)
```

## Arguments

- n_samples:

  Number of samples (donors).

- cells_per_sample:

  Number of cells simulated per sample.

- n_genes:

  Number of genes.

- n_celltypes:

  Number of cell types (named `CT1`, `CT2`, ...).

- prop_responsive:

  Proportion of genes with a true effect for each exposure.

- effect_sd:

  Standard deviation of the true (log-scale) exposure effects.

- prop_celltype_responsive:

  Probability that a cell type responds to a given exposure. `CT1`
  always responds to every exposure.

- size:

  Negative binomial size (dispersion) parameter. Smaller values give
  more overdispersion.

## Value

A `SingleCellExperiment` with:

- `assays$counts`:

  Genes x cells count matrix.

- `colData`:

  `sample`, `pm25`, `bpa`, `smoker`, `age`, `celltype`.

- `rowData`:

  `beta_pm25`, `beta_bpa`, `beta_smoker`: true effects, where nonzero
  means a truly responsive gene.

- `metadata$responsive`:

  Cell type x exposure matrix of which cell types respond to which
  exposures.

## Details

Exposures are constant within a sample, so the effective sample size is
the number of samples, not the number of cells.

## Examples

``` r
sce <- make_example_sc()
```
