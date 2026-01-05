# tidyexposomics 0.99.11

* Setting `nipalsMCIA` to suggests due to `nipalsMCIA` RMD check failures.

# tidyexposomics 0.99.10

* Fixed `set.seed()` to make network plots consistent within vignette.
* Fixed issue with impute exposure logic.
* Updated naming conventions for variable maps.

# tidyexposomics 0.99.9

* Added `set.seed()` to make network plots consistent.
* Updated vignette formatting and PCA documentation.

# tidyexposomics 0.99.8

* Added BiocFileCache-based caching and removed globalenv assignments.
* Swithced to small example data for vignette and added figure captions.
* Added MultiAssayExperiment to Depends.
* Standardized MIT license file.
* Updated pipeline summary output.

# tidyexposomics 0.99.5/0.99.6/0.99.7

* Removed .Rproj tracking, .gitignore lines.
* Updated testing strategy.
* Modified dependency to 4.5.0.
* Added vignette chunk labels.

# tidyexposomics 0.99.2/0.99.3/0.99.4

* Updated namespace issues with tidybulk.

# tidyexposomics 0.99.1

* Updated variable language for consistency.
* Pre-build objects for vignette run time.
* Updated circos plot visuals.
* Updated network functions to create unique node names.


# tidyexposomics 0.99.0

* Initial Bioconductor submission.
* Implements a full exposome-omics pipeline including:
  - Quality control 
  - Association analysis
  - Differential abundance analysis
  - Multi-omics integration analysis
  - Functional enrichment analysis
  - Vizualization 


