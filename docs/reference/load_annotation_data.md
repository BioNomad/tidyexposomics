# Load Ontology Data

Downloads and loads `chebi`, `ecto`, and `hpo`.

## Usage

``` r
load_annotation_data(
  to_load = c("all", "chebi", "ecto", "hpo"),
  verbose = TRUE,
  validate = TRUE
)
```

## Arguments

- to_load:

  Character vector indicating which ontology to load.

- verbose:

  Logical; print messages.

- validate:

  Logical; validate MD5 checksum.

## Value

A named list of ontology objects.

## Examples

``` r
if (FALSE) { # \dontrun{
# load single ontology
onts <- load_annotation_data(to_load = "ecto")

# load all ontologies
onts <- load_annotation_data(to_load = "all")
} # }
```
