# tidyexposomics <a href="#"><img src="./inst/logo.png" align="right" height="200" /></a>

<br>
<br>
<br>
<br>

## Overview

The `tidyexposomics` package is designed to facilitate the integration of exposure and omics data to identify exposure-omics associations. We structure our commands to fit into the tidyverse framework, where commands are designed to be simplified and intuitive. Here we provide functionality to perform quality control, sample and exposure association analysis, differential abundance analysis, multi-omics integration, and functional enrichment analysis.

<p align="center" width="100%">
    <img width="80%" src="./inst/overview.png"> 
</p>

## Installation

You can install the development version of `tidyexposomics` from GitHub with:

```r
# Install directly through GitHub
devtools::install_github("BioNomad/tidyexposomics")
```

## Command Structure

To make the package more user-friendly, we have named our functions to be more intuitive. For example, we use the following naming conventions:

- `associate_*` Associations using models
- `correlate_*` Correlation based approaches
- `run_*` Running more complex pipelines
- `plot_*` Visualizing results
- `filter_*`, `transform_*`, `extract_*`, and `pivot_*` Data pre-processing steps

We provide functionality to either add results to the existing object storing the omics/exposure data or to return the direct results using the get option. We suggest adding results, as we also include a run_bibliography function to generate a bibliography of the workflow. This is useful for keeping track of the pipeline steps.

## Tutorial

`tidyexposomics` tutorials and documentation can be found at:

[LINK TO BE INSERTED LATER]
