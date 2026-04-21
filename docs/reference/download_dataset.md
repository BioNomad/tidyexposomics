# Download and cache a tidyexposomics dataset

Download and cache a tidyexposomics dataset

## Usage

``` r
download_dataset(
  name = c("omics_list", "fdata", "meta", "annotated_cb", "expom_1", "chebi", "ecto",
    "hpo"),
  verbose = TRUE,
  validate = TRUE
)
```

## Arguments

- name:

  Dataset name: one of "omics_list", "fdata", "meta", "annotated_cb",
  "expom_1", "chebi", "ecto", "hpo".

- verbose:

  Logical; print messages.

- validate:

  Logical; validate MD5 checksum.

## Value

A list or object loaded from the cached .RData.
