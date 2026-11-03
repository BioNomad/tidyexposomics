# Make an example codebook for the simulated single-cell data

Returns a data dictionary describing the `colData()` variables produced
by
[`make_example_sc()`](https://BioNomad.github.io/tidyexposomics/reference/make_example_sc.md).
Each row is one variable, annotated with a label, its role in an
analysis, its exposure domain, the level at which it varies, and its
unit.

## Usage

``` r
make_example_sce_codebook()
```

## Value

A data frame with one row per variable and columns:

- `variable`:

  Column name in `colData()`.

- `label`:

  Human-readable description.

- `role`:

  One of `"id"`, `"exposure"`, `"covariate"`, or `"cell_annotation"`.

- `domain`:

  Category of the variable (e.g. `"air pollution"`).

- `level`:

  Whether the variable varies by `"sample"` or by `"cell"`.

- `unit`:

  Unit of measurement, or `NA` if not applicable. Units are
  illustrative, since the simulator does not assign real ones.

## Examples

``` r
make_example_sce_codebook()
#>   variable                           label            role        domain  level
#> 1   sample                       Sample ID              id    identifier sample
#> 2     pm25 Fine particulate matter (PM2.5)        exposure air pollution sample
#> 3      bpa               Bisphenol A (BPA)        exposure      chemical sample
#> 4   smoker          Current smoking status        exposure     lifestyle sample
#> 5      age               Age at collection       covariate   demographic sample
#> 6 celltype            Cell type annotation cell_annotation     cell type   cell
#>    unit
#> 1  <NA>
#> 2 ug/m3
#> 3 ng/mL
#> 4  <NA>
#> 5 years
#> 6  <NA>

make_example_sce_codebook() |>
    dplyr::filter(role == "exposure")
#>   variable                           label     role        domain  level  unit
#> 1     pm25 Fine particulate matter (PM2.5) exposure air pollution sample ug/m3
#> 2      bpa               Bisphenol A (BPA) exposure      chemical sample ng/mL
#> 3   smoker          Current smoking status exposure     lifestyle sample  <NA>
```
