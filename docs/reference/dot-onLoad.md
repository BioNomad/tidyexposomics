# Internal - onLoad hook to register www assets

Registers `inst/app/www` as a Shiny resource path `"www"` if present.

## Usage

``` r
.onLoad(libname, pkgname)
```

## Value

Invisibly returns `NULL`
