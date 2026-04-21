# Build the Ontology Annotation Shiny app

Returns a Shiny app object for ontology annotation. This function does
not launch the app. Please call
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) yourself
(see examples).

## Usage

``` r
build_ont_annot_app(use_demo = TRUE, ...)
```

## Arguments

- use_demo:

  Logical, if TRUE, load packaged lightweight demo ontologies.

- ...:

  Optional named overrides passed as data.frames with columns `id`,
  `name`, `ancestors`, `parents`, `children`: `hpo`, `ecto`, `chebi`.

## Value

A `shiny.appobj`.

## Examples

``` r
if (interactive()) {
    app <- build_ont_annot_app()
    shiny::runApp(app)
}
```
