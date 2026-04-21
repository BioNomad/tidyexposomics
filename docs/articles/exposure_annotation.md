# Exposure Annotation

## Exposure Metadata and Ontology Annotation

------------------------------------------------------------------------

### Codebook Setup

Before starting an exposomics data analysis we recommend having a
codebook, with information on your exposure variables. Some suggestions:

- **Variable Name**: The name of the variable in the data set.

- **Variable Description**: A concise description of what the variable
  measures, including units (e.g., “urinary bisphenol A (ng/mL)”).

- **Variable Type**: The type of variable, such as continuous,
  categorical, or binary.

- **Variable Period**: The period of time over which the variable was
  measured, such as “lifetime”, “year”, “month”, or “day”.

- **Variable Location**: The location where the variable was measured,
  such as “home”, “work”, “school”, or “geospatial code”.

- **Variable Ontology**: The ontology term associated with the variable.

### Ontology Choices

Variables captured in the codebook should be annotated with ontology
terms to provide a standardized vocabulary for the variables. We
recommend using the following ontologies for exposure and outcome
variables:

- [Environment Exposure
  Ontology](https://www.ebi.ac.uk/ols4/ontologies/ecto) to annotate your
  exposure variables.

- [Human Phenotype Ontology](https://www.ebi.ac.uk/ols4/ontologies/hp)
  to annotate your outcome variables and phenotypic data.

- [Chemical Entities of Biological
  Interest](https://www.ebi.ac.uk/ols4/ontologies/chebi) to annotate
  your chemical exposure variables.

**Why annotate with ontologies?**

- **Interpretability**: Ontology labels clarify ambiguous or
  inconsistently named variables.

- **Harmonization**: You can compare and combine variables across
  datasets when they map to the same term.

- **Grouping**: Ontologies allow you to collapse fine-grained exposures
  into broader categories.

- **Integration**: Many public tools, knowledge graphs, and repositories
  are ontology-aware. This can make your results more interoperable and
  reusable.

### Ontology Annotation App

To help annotate exposure variables, we provide a lightweight shiny app:

``` r
# Launch the shiny app to annotate exposure variables
shiny::runApp(build_ont_annot_app())
```

**To use the app:**

- Click **Browse** to select your exposure metadata file.

- Then you can click the variable you’d like to link to an annotation
  term and search in the **Choose Ontology Term** dropdown.

- After you select a term, you will see a short description of the term.

- After you are done, click **Apply Annotate** to save the annotation.

- Now you can group exposures into larger categories by selecting each
  line and then choose your ontology and root depth level (where a lower
  number means a more general term).

- Then you can click **Apply Categorization** to apply the selected
  categorization to the selected rows.

- If the ontology has nothing to do with your variable, you may manually
  enter a category in the **Category** column. This will change the
  **Category Source** to manual and will not be linked to the ontology.

- Once you have annotated all your variables click **Download Annotated
  CSV** to save the annotated metadata file.

![Screenshot of the ontology annotation app. The sidebar has an upload
button to load your exposure metadata file. The main panel displays the
uploaded exposure metadata where users can select variables to annotate
and categorize. The sidebar also contains buttons to apply
annotations/categorizations and download the annotated
file.](ont_annot.png)

Screenshot of the ontology annotation app. The sidebar has an upload
button to load your exposure metadata file. The main panel displays the
uploaded exposure metadata where users can select variables to annotate
and categorize. The sidebar also contains buttons to apply
annotations/categorizations and download the annotated file.

### Session Info

See Session Info

``` r

sessionInfo()
## R version 4.5.1 (2025-06-13 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26200)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] BiocStyle_2.38.0
## 
## loaded via a namespace (and not attached):
##  [1] cli_3.6.5           knitr_1.51          rlang_1.1.7        
##  [4] xfun_0.54           otel_0.2.0          textshaping_1.0.4  
##  [7] jsonlite_2.0.0      htmltools_0.5.9     ragg_1.5.0         
## [10] sass_0.4.10         rmarkdown_2.30      evaluate_1.0.5     
## [13] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
## [16] lifecycle_1.0.5     bookdown_0.46       BiocManager_1.30.27
## [19] compiler_4.5.1      fs_2.1.0            htmlwidgets_1.6.4  
## [22] rstudioapi_0.18.0   systemfonts_1.3.2   digest_0.6.39      
## [25] R6_2.6.1            bslib_0.10.0        tools_4.5.1        
## [28] pkgdown_2.2.0       cachem_1.1.0        desc_1.4.3
```
