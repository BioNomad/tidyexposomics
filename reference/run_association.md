# Run Association Analysis

Perform GLM-based association testing between a specified outcome and
features from exposures, omics, or latent factors. Each feature is
tested in its own model, built from a user-supplied formula. Models can
be run on subsets of samples, stratified by a grouping variable, or
repeated across several formulas. Every fitted model is stored under a
name so that specifications can be compared later with
[`compare_association()`](https://BioNomad.github.io/tidyexposomics/reference/compare_association.md).

## Usage

``` r
run_association(
  exposomicset,
  formula,
  source = c("omics", "exposures", "factors"),
  feature_set = NULL,
  log_trans = TRUE,
  top_n = NULL,
  family = "gaussian",
  correction_method = "fdr",
  subset = NULL,
  strata = NULL,
  name = NULL,
  overwrite = FALSE,
  action = "add",
  feature_col = NULL,
  mirna_assays = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing data and metadata.

- formula:

  A two-sided formula whose left-hand side is a single column of
  `colData` (the outcome) and whose right-hand side contains the
  covariates, e.g. `smoker ~ age * sex`. Each feature is added to the
  model as a main effect. To control how the feature enters the model
  (for example to test an interaction), use the placeholder `.feature`
  on the right-hand side, e.g. `smoker ~ .feature * sex + age`. Use
  `outcome ~ 1` for an unadjusted model. Alternatively, a **named list
  of formulas** (e.g. `list(crude = y ~ 1, adjusted = y ~ age + sex)`)
  fits every specification in one call, each stored under its list name.

- source:

  Source of features to test. One of `"omics"`, `"exposures"`,
  `"factors"`.

- feature_set:

  Optional character vector of exposure or GO terms to test.

- log_trans:

  Optional boolean value dictating whether or not to log2 transform
  omics features.

- top_n:

  Optional integer: if using omics source, select top `n` most variable
  features.

- family:

  GLM family; `"gaussian"` or `"binomial"`.

- correction_method:

  Method for p-value adjustment (default: `"fdr"`). Adjustment is
  performed separately for each model and each effect (see Value).

- subset:

  Optional bare expression evaluated against the (unscaled) `colData`
  that selects the samples to model, e.g. `sex == "female"` or
  `age > 40 & smoker == "no"`.

- strata:

  Optional single column name of `colData`. The model is fit separately
  within each level and stored as `[MODEL_NAME]_[LEVEL_IN_STRATA]`. The
  stratifying variable must not appear in `formula`.

- name:

  Optional name for the stored model. Defaults to `"default"`, or to the
  text of `subset` when one is given. Ignored when `formula` is a named
  list (the list names are used).

- overwrite:

  Logical; if `TRUE`, replace an existing model with the same name
  instead of erroring or appending features to it (default: `FALSE`).

- action:

  If `"add"` (default), saves results to metadata; else returns results
  as list.

- feature_col:

  The column in `rowData` for matching gene symbols or IDs.

- mirna_assays:

  Optional character vector of assays to exclude when extracting GO
  terms.

## Value

If `action = "add"`, returns updated `MultiAssayExperiment`. Otherwise,
returns a list of:

- `models`: a named list of fitted models, each containing `results_df`
  (tidy summary of associations), `formula`, `covariates`, `family`,
  `subset`, `strata_var`, `stratum`, and `n_samples`

- `model_data`: the (full) data used in the GLMs

In `results_df`, `term` is the feature, `coef` is the model coefficient
name, `effect` is `"main"` for the feature's main effect or the name of
the interacting term (e.g. `"sexMale"`), and `n` / `n_events` give the
observations (and binomial events) used.

## Details

Numeric `colData` columns are z-scored once on the **full** sample
before any subsetting, so effect sizes are on the same standardized
scale in every subset and stratum and can be compared directly.

If a model with the same `name` already exists for the source and was
fit with the same formula, subset, and family, new features are appended
to it and p-values are re-adjusted across the combined set. A different
specification under an existing name is an error unless
`overwrite = TRUE`.

With `strata`, a stratum that cannot be fit (for example, a single
outcome level or a constant covariate) is skipped with a warning.

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

# covariate-adjusted main effect
mae <- mae |>
    run_association(
        formula = smoker ~ age,
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial"
    )
#> Running GLMs.

# several specifications at once
mae <- mae |>
    run_association(
        formula = list(
            crude = smoker ~ 1,
            adjusted = smoker ~ age
        ),
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial"
    )
#> Running GLMs.

if (FALSE) { # \dontrun{
# subset of samples (assumes a `sex` column in colData)
mae <- mae |>
    run_association(
        formula = smoker ~ age,
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial",
        subset = sex == "female",
        name = "female"
    )

# stratified by sex, then compared
mae <- mae |>
    run_association(
        formula = smoker ~ age,
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        family = "binomial",
        strata = "sex",
        name = "by_sex"
    )

mae |>
    compare_association(
        source = "exposures",
        models = c("by_sex_female", "by_sex_male")
    )
} # }
```
