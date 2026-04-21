# Internal - categorizes to ontology roots / depth

Internal - categorizes to ontology roots / depth

## Usage

``` r
.run_categorize_ontology(
  data,
  id_col,
  ontologyDF,
  root_level = 0,
  assign_label = TRUE
)
```

## Arguments

- data:

  Data frame with an ID column.

- id_col:

  Column name containing ontology term IDs.

- ontologyDF:

  Processed ontology data.frame (with list columns for relationships and
  `depth`).

- root_level:

  Either a character vector of explicit root IDs, a numeric depth, or
  (default) top-level roots.

- assign_label:

  If TRUE, return input with new columns; otherwise return assigned IDs.

## Value

Data frame with `root_id`, `root_label`, `category`, `category_source`
(if `assign_label`).
