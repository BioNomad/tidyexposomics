# Introduction to tidyexposomics

Exposure to environmental factors is a major determinant of health and
disease. The exposome is a term that represents the totality of
environmental exposures (or in practical terms, the integrated
compilation of all physical, chemical, biological, and (psycho)social
influences that impact biology) from conception to death, that drive
disease and overall health (Miller & Banbury Exposomics Consortium,
2025; Wild, 2005). Characterizing the relationship between a
multiplicity of exposures and disease outcomes can be a daunting task.
Epidemiological studies have moved towards exposure-wide association
studies (ExWAS) where sets of exposures are associated with an outcome
(Chung et al., 2024). The advent of high-throughput technologies has
enabled profiling of several layers of biological information. These
layers can also be integrated to better understand the relationship
between exposures and disease outcomes.

The `tidyexposomics` package is designed to facilitate the integration
of exposure and omics data to identify exposure-omics associations and
their relevance to health outcomes.`tidyexposomics` extends the
tidy-Bioconductor ecosystem (e.g., tidybulk, tidySummarizedExperiment)
to exposome multi-omics integration using the MultiAssayExperiment
container. It provides tidyverse-style accessors and functions for
association testing, multi-omics integration, and ontology-driven
enrichment, in an effort to complement existing tidy-Bioc tools.

![tidyexposomics pipeline overview. QC, association testing,
integration, and enrichment steps on a
MultiAssayExperiment.](overview.png)

tidyexposomics pipeline overview. QC, association testing, integration,
and enrichment steps on a MultiAssayExperiment.

## Installation

``` r
# install the package
BiocManager::install("tidyexposomics")

# load the package
library(tidyexposomics)
```

## Command Structure

------------------------------------------------------------------------

To make the package more user-friendly, we have named our functions to
be more intuitive. For example, we use the following naming conventions:

![Command naming conventions used throughout \`tidyexposomics\`. More
complex pipelines begin with the \`run\_\` prefix, visualizations with
\`plot\_\`, and data processing with \`filter\_\`, \`transform\_\`,
\`pivot\_\`, or \`extract\_\` prefixes.](command_str.png)

Command naming conventions used throughout `tidyexposomics`. More
complex pipelines begin with the `run_` prefix, visualizations with
`plot_`, and data processing with `filter_`, `transform_`, `pivot_`, or
`extract_` prefixes.

We provide functionality to either add results to the existing object
storing the omics/exposure data or to return results directly using
`action = "get"`. We suggest adding results, given that pipeline steps
are tracked and can be output to the R console, plotted as a workflow
diagram, or exported to an Excel worksheet.

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

## Loading Data

------------------------------------------------------------------------

To get started we need to load the data. The `create_exposomicset`
function is used to create a `MultiAssayExperiment` object that contains
exposure and omics data. As a quick introduction, a
`MultiAssayExperiment` object is a container for storing multiple assays
(e.g., omics data) and their associated metadata:

![Overview of the MultiAssayExperiment structure linking samples,
assays, and metadata.](mae_rep.png)

Overview of the MultiAssayExperiment structure linking samples, assays,
and metadata.

We use the MultiAssayExperiment object to store the exposure and omics
data. The create_expomicset function has several arguments:

- `codebook`: is a data frame that contains information about the
  variables in the exposure metadata. The column names must contain
  variable where the values are the column names of the exposure data
  frame, and category which contains general categories for the variable
  names. This is the data frame you created with the ontology annotation
  app!

- `exposure`: is a data frame that contains the exposure and other
  metadata.

- `omics`: is a list of data frames that contain the omics data.

- `row_data`: argument is a list of data frames that contain information
  about the rows of each omics data frame.

We are going to start by loading in example data pulled from the
[ISGlobal Exposome Data Challenge
2021](https://doi.org/10.1016/j.envint.2022.107422) (Maitre et al.,
2022). Specifically, we will examine how exposures and omics features
relate to asthma status in asthma patients with a lower socioeconomic
status (SES).

``` r
# Load Libraries
library(tidyverse)
library(tidyexposomics)
```

``` r
# Load example data
data("tidyexposomics_example")

# Create exposomic set object
expom <- create_exposomicset(
    codebook = tidyexposomics_example$annotated_cb,
    exposure = tidyexposomics_example$meta,
    omics = list(
        "Gene Expression" = tidyexposomics_example$exp_filt,
        "Methylation" = tidyexposomics_example$methyl_filt
    ),
    row_data = list(
        "Gene Expression" = tidyexposomics_example$exp_fdata,
        "Methylation" = tidyexposomics_example$methyl_fdata
    )
)
```

    ## Ensuring all omics datasets are matrices with column names.

    ## Creating SummarizedExperiment objects.

    ## Creating MultiAssayExperiment object.

    ## MultiAssayExperiment created successfully.

We are interested in how the exposome affects health outcomes, so let’s
define which metadata variables represent exposure variables.

``` r
# Grab exposure variables
exp_vars <- tidyexposomics_example$annotated_cb |>
    filter(category %in% c(
        "aerosol",
        "main group molecular entity",
        "polyatomic entity"
    )) |>
    pull(variable) |>
    as.character()
```

## Quality Control

------------------------------------------------------------------------

### Missingness

Oftentimes when collecting data, there are missing values. Let’s use the
`plot_missing` function to determine where our missing values are:

``` r
# Plot the missingness summary
plot_missing(
    exposomicset = expom,
    plot_type = "summary",
    threshold = 0
)
```

![Count of features with missing data above a 0% missingness threshold
by data layer. Exposure data have variables with
missingness.](tidyexposomics_files/figure-html/missing-bar-1.png)

Count of features with missing data above a 0% missingness threshold by
data layer. Exposure data have variables with missingness.

Here we see that there are 4 variables in the exposure data that are
missing data. Let’s take a look at them:

``` r
# Plot missing variables withing exposure group
plot_missing(
    exposomicset = expom,
    plot_type = "lollipop",
    threshold = 0,
    layers = "Exposure"
)
```

![Percent missingness per exposure variable. Parity, \`h_parity_None\`,
shows the highest
missingness.](tidyexposomics_files/figure-html/missing-bar-lollipop-1.png)

Percent missingness per exposure variable. Parity, `h_parity_None`,
shows the highest missingness.

Here we see that one variable, `h_parity_None`, has about 4% missing
values. We can apply a missingness filter using the `filter_missing`
function. However, given that this level of missingness is quite low, we
will not be applying a missingness filter and instead impute the missing
data.

### Imputation

The `run_impute_missing` function is used to impute missing values. Here
we can specify the imputation method for exposure and omics data
separately.

The `exposure_impute_method` argument is used to set the imputation
method for exposure data, and the `omics_impute_method` argument is used
to set the imputation method for omics data. The `omics_to_impute`
argument is used to specify which omics data to impute. Here we will
impute the exposure data given using the `missforest` method, but other
options for imputation methods include:

- `median`: Imputes missing values with the median of the variable.

- `mean`: Imputes missing values with the mean of the variable.

- `knn`: Uses k-nearest neighbors to impute missing values.

- `mice`: Uses the Multivariate Imputation by Chained Equations (MICE)
  method to impute missing values.

- `missforest`: Uses the MissForest method to impute missing values.

- `lod_sqrt2`: Imputes missing values using the square root of the lower
  limit of detection (LOD) for each variable. This is useful for
  variables that have a lower limit of detection, such as chemical
  exposures.

``` r
# Impute missing values
expom <- run_impute_missing(
    exposomicset = expom,
    exposure_impute_method = "missforest",
    exposure_cols = exp_vars
)
```

    ## Imputing exposure data using method: missforest

### Filtering Omics Features

We can filter omics features based on variance or expression levels. The
`filter_omics` function is used to filter omics features. The method
argument is used to set the method for filtering. Here we can use
either:

- **Variance**: Filters features based on variance. We recommend this
  for omics based on continuous measurements, such as log-transformed
  counts, M-values, protein intensities, or metabolite concentrations.

- **Expression**: Filters features based on expression levels. We
  recommend this for omics where many values may be near-zero or zero,
  such as RNA-seq data.

The `assays` argument is used to specify which omics data to filter. The
`assay_name` argument is used to specify which assay to filter. The
`min_var`, `min_value`, and `min_prop` arguments are used to set the
minimum variance, minimum expression value, and minimum proportion of
samples exceeding the minimum value, respectively.

``` r
# filter omics layers by variance and expression
# Methylation filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "variance",
    assays = "Methylation",
    assay_name = 1,
    min_var = 0.05
)
```

    ## Filtering assay: Methylation

    ## Filtered 223 of 500 features from 'Methylation' using method 'variance'

``` r
# Gene expression filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "expression",
    assays = "Gene Expression",
    assay_name = 1,
    min_value = 1,
    min_prop = 0.3
)
```

    ## Filtering assay: Gene Expression

    ## Filtered 29 of 500 features from 'Gene Expression' using method 'expression'

### Normality Check

When determining variable associations, it is important to check the
normality of the data. The `run_normality_check` function is used to
check the normality of the data.

``` r
# Check variable normality
expom <- run_normality_check(
    exposomicset = expom,
    action = "add"
)
```

    ## Checking Normality Using Shapiro-Wilk Test

    ## 9 Exposure Variables are Normally Distributed

    ## 6 Exposure Variables are NOT Normally Distributed

The `transform_exposure` function is used to transform the data to make
it more normal. Here the transform_method is set to `boxcox_best` as it
will automatically select the best transformation method based on the
data. The `transform_method` can be manually set to `log2`, `sqrt`, or
`x_1_3` as well. We specify the exposure_cols argument to set the
columns to transform.

``` r
# Transform variables
expom <- transform_exposure(
    exposomicset = expom,
    transform_method = "boxcox_best",
    exposure_cols = exp_vars
)
```

    ## Applying the boxcox_best transformation.

To check the normality of the exposure data, we can use the
`plot_normality_summary` function. This function plots the normality of
the data before and after transformation. The `transformed` argument is
set to `TRUE` to plot the normality status of the transformed data.

``` r
# Examine normality summary
plot_normality_summary(
    exposomicset = expom,
    transformed = TRUE
)
```

![Normality status of numeric exposure variables after Box-Cox
transformation.](tidyexposomics_files/figure-html/norm-plot-1.png)

Normality status of numeric exposure variables after Box-Cox
transformation.

### Principal Component Analysis

To identify the variability of the data, we can perform a principal
component analysis (PCA). The `run_pca` performs a joint PCA across all
numeric exposures and omic assays after standardization, identifying
shared axes of variation across layers. The resulting PCs in `colData()`
reflect integrated sample-level variance across all data types, and
outliers are defined in that joint multi-omics PC space.

Here we specify that we would like to log-transform the exposure and
omics data before performing PCA using the `log_trans_exp` and the
`log_trans_omics` arguments, respectively. We automatically identify
sample outliers based on the Mahalanobis distance, a measure of the
distance between a point and a distribution.

``` r
# Perform principal component analysis
expom <- run_pca(
    exposomicset = expom,
    log_trans_exp = TRUE,
    log_trans_omics = TRUE,
    action = "add"
)
```

    ## Identifying common samples.

    ## Subsetting exposure data.

    ## Subsetting omics data.

    ## Performing PCA on Feature Space.

    ## Performing PCA on Sample Space.

    ## Outliers detected: s1231

``` r
# Plot the PCA plot of sample and feature space
plot_pca(exposomicset = expom)
```

![PCA of sample and feature space with sample outlier
detection.](tidyexposomics_files/figure-html/pca-plot-1.png)

PCA of sample and feature space with sample outlier detection.

Here we see one sample outlier, and that most variation is captured in
the first two principal components for both features and samples. We can
filter out the outlier using the `filter_sample_outliers` function.

``` r
# Filter out sample outliers
expom <- filter_sample_outliers(
    exposomicset = expom,
    outliers = c("s1231")
)
```

    ## Removing outliers: s1231

To understand the relationship between the principal components and
exposures we can correlate them using the `run_correlation` function.
Here we specify that the `feature_type` is `pcs` for principal
components, specify a set of exposure variables, `exp_vars`, and the
number of principal components, `n_pcs.` We set `correlation_cutoff` to
`0` and `pval_cutoff` to `1` to initially include all correlations.

``` r
# Run the correlation analysis
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "pcs",
    exposure_cols = exp_vars,
    n_pcs = 20,
    action = "add",
    correlation_cutoff = 0,
    pval_cutoff = 1
)
```

We can visualize these correlations with the `plot_correlation_tile`
function. We specify we are plotting the `feature_type` of `pcs` to grab
the principal component correlation results. We then set the
significance threshold to 0.05 with the `pval_cutoff` argument.

``` r
# Plot the correlation tile plot
plot_correlation_tile(
    exposomicset = expom,
    feature_type = "pcs",
    pval_cutoff = 0.05
)
```

![Correlation heatmap of exposures versus principal components. Child
lead levels (\`hs_pb_c_Log2\`) and maternal BPA levels
(\`hs_bpa_madj_Log2\`) are associated with the most principal
components.](tidyexposomics_files/figure-html/plot-pc-corr-1.png)

Correlation heatmap of exposures versus principal components. Child lead
levels (`hs_pb_c_Log2`) and maternal BPA levels (`hs_bpa_madj_Log2`) are
associated with the most principal components.

### Exposure Summary

We can summarize the exposure data using the `run_summarize_exposures`
function. This function calculates summary statistics for each exposure
variable, including the number of values, number of missing values,
minimum, maximum, range, sum, median, mean, standard error, and
confidence intervals. The `exposure_cols` argument determines which
variables to include in the summary.

``` r
# Summarize exposure data
run_summarize_exposures(
    exposomicset = expom,
    action = "get",
    exposure_cols = exp_vars
) |>
    head()
```

    ## # A tibble: 6 × 27
    ##   variable   n_values  n_na   min   max range    sum median  mean    se ci_lower
    ##   <chr>         <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>    <dbl>
    ## 1 h_pm10_ra…       47     0 16.6  25.8   9.25 989.    21.2  21.0   0.32    20.4 
    ## 2 h_pm25_ra…       47     0 10.6  18.2   7.55 697.    14.6  14.8   0.25    14.3 
    ## 3 hs_bpa_ma…       47     0  0     2.31  2.31  71.3    1.58  1.52  0.06     1.41
    ## 4 hs_mibp_c…       47     0  0.1   0.22  0.11   7.65   0.17  0.16  0        0.16
    ## 5 hs_pb_c_L…       47     0  1.44  5     3.56 156.     3.31  3.31  0.11     3.1 
    ## 6 hs_pfhxs_…       47     0  0     2.14  2.13  65.7    1.48  1.4   0.07     1.25
    ## # ℹ 16 more variables: ci_upper <dbl>, variance <dbl>, sd <dbl>,
    ## #   coef_var <dbl>, period <chr>, location <chr>, description <chr>,
    ## #   var_type <chr>, transformation <chr>, selected_ontology_label <chr>,
    ## #   selected_ontology_id <chr>, root_id <chr>, root_label <chr>,
    ## #   category <chr>, category_source <chr>, transformation_applied <chr>

### Exposure Visualization

------------------------------------------------------------------------

To visualize our exposure data, we can use the `plot_exposures`
function. This function allows us to plot the exposure data in a variety
of ways. Here we will plot the exposure data using a boxplot. The
`exposure_cat` argument is used to set the exposure category to plot.
Additionally, we could specify `exposure_cols` to only plot certain
exposures. The `plot_type` argument is used to set the type of plot to
create. Here we use a boxplot, but we could also use a ridge plot.

``` r
# Plot aerosol exposure distributions by sex
plot_exposures(
    exposomicset = expom,
    group_by = "e3_sex_None",
    exposure_cat = "aerosol",
    plot_type = "boxplot",
    ylab = "Values",
    title = "Aerosol Exposure by Sex"
)
```

![Distribution of aerosol exposures by
sex.](tidyexposomics_files/figure-html/plot-aerosol-1.png)

Distribution of aerosol exposures by sex.

Here we do not see any significant differences in aerosol exposure
between males and females.

## Sample-Exposure Association

------------------------------------------------------------------------

### Sample Clustering

The `run_cluster_samples` function is used to cluster samples based on
the exposure data, clustering approaches are available by setting the
`clustering_approach` argument. Here we use the `dynamic` approach,
which uses a dynamic tree cut method to identify clusters. Other options
are:

- `gap`: **Gap statistic method (default)**; estimates optimal `k` by
  comparing within-cluster dispersion to that of reference data.

- `diana`: **Divisive hierarchical clustering (DIANA)**; chooses `k`
  based on the largest drop in dendrogram height.

- `elbow`: **Elbow method**; detects the point of maximum curvature in
  within-cluster sum of squares (WSS) to determine `k`.

- `dynamic`: **Dynamic tree cut**; adaptively detects clusters from a
  dendrogram structure without needing to predefine `k`.

- `density`: **Density-based clustering** (via `densityClust`);
  identifies clusters based on local density peaks in distance space.

``` r
# Sample clustering
expom <- run_cluster_samples(
    exposomicset = expom,
    exposure_cols = exp_vars,
    clustering_approach = "dynamic",
    action = "add"
)
```

    ## Starting clustering analysis...

    ##  ..cutHeight not given, setting it to 40.7  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

    ## Optimal number of clusters for samples: 2

We plot the sample clusters using the `plot_sample_clusters` function.
This function plots z-scored values of the exposure data for each
sample, colored by the cluster assignment. The `exposure_cols` argument
is used to set the columns to plot.

``` r
# Plot the sample clusters
plot_sample_clusters(
    exposomicset = expom,
    exposure_cols = exp_vars
)
```

    ## tidyHeatmap says: If you use tidyHeatmap for scientific research, please cite: Mangiola, S. and Papenfuss, A.T., 2020. 'tidyHeatmap: an R package for modular heatmap production based on tidy principles.' Journal of Open Source Software. doi:10.21105/joss.02472.
    ## This message is displayed once per session.

    ## Warning: `when()` was deprecated in purrr 1.0.0.
    ## ℹ Please use `if` instead.
    ## ℹ The deprecated feature was likely used in the tidyHeatmap package.
    ##   Please report the issue at
    ##   <https://github.com/stemangiola/tidyHeatmap/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![Sample clustering heatmap using exposure profiles (z-scored). Clusters
appear mostly driven by aerosol exposure during
pregnancy.](tidyexposomics_files/figure-html/plot-sample-clusters-1.png)

Sample clustering heatmap using exposure profiles (z-scored). Clusters
appear mostly driven by aerosol exposure during pregnancy.

Here we see two clusters, largely driven by particulate matter/aerosol
exposure during pregnancy (`h_pm25_ratio_preg_None` and
`h_pm10_ratio_preg_None`).

### Exposure Correlations

The `run_correlation` function identifies correlations between exposure
variables. We set `feature_type` to `exposures` to focus on exposure
variables and use a correlation cutoff of `0.3` to filter for meaningful
associations. This cutoff can be adjusted based on your data and
analysis needs.

``` r
# Run correlation analysis
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "exposures",
    action = "add",
    exposure_cols = exp_vars,
    correlation_cutoff = 0.3
)
```

To visualize the exposure correlations, we can use the
`plot_circos_correlation` function. Here we will plot the circos plot.
This function creates a circular plot of the exposure correlations. The
`correlation_cutoff` argument is used to set the minimum correlation
score for the association. Here we use a cutoff of `0.3`.

``` r
# Plot exposure correlation circos plot
plot_circos_correlation(
    exposomicset = expom,
    feature_type = "exposures",
    corr_threshold = 0.3,
    exposure_cols = exp_vars
)
```

![Circos view of exposure-exposure correlations (threshold
0.3).](tidyexposomics_files/figure-html/exposure-circos-corr-1.png)

Circos view of exposure-exposure correlations (threshold 0.3).

### Exposure-wide association (ExWAS)

The `run_association` function performs an ExWAS analysis to identify
associations between exposures and outcomes. We specify the data source,
outcome variable, feature set, and covariates for the analysis. Since we
have a binary outcome, we set the model family to binomial.

``` r
# Perform ExWAS
expom <- run_association(
    exposomicset = expom,
    source = "exposures",
    outcome = "hs_asthma",
    feature_set = exp_vars,
    action = "add",
    family = "binomial"
)
```

    ## Running GLMs.

To visualize the results of the ExWAS analysis, we can use the
`plot_association` function, which will plot results for the the
specified features. The `terms` argument is used to set the features to
plot. The `filter_thresh` argument is used to set the threshold for
filtering the results. The `filter_col` argument is used to set the
column to filter on. Here we use `p.value` to filter on the p-value of
the association. We can also include the R^2 or adjusted R^2 (if
covariates are included) using the `r2_col` argument.

``` r
# Plot the association forest plot
plot_association(
    exposomicset = expom,
    source = "exposures",
    terms = exp_vars,
    filter_thresh = 0.05,
    filter_col = "p.value",
    r2_col = "r2"
)
```

![ExWAS associations of exposures with asthma status. No exposures are
significantly associated (P \< 0.05) with asthma
status.](tidyexposomics_files/figure-html/plot-exposure-assoc-1.png)

ExWAS associations of exposures with asthma status. No exposures are
significantly associated (P \< 0.05) with asthma status.

Here we see that no exposure variables are significantly associated with
our asthma status. Although we do see that confidence interval for child
Mono-iso-butyl phthalate (MiBP) levels (`hs_mibp_cadj_Log2`) does not
cross 0, indicating a negative, albeit not significant (P \< 0.05)
association.

We can also associate our omics features with an outcome of interest
using the `run_association` function. Here we specify an additional
argument, `top_n`, which is used to set the top number of high variance
omics features to include per omics layer.

``` r
# Perform ExWAS
expom <- run_association(
    exposomicset = expom,
    outcome = "hs_asthma",
    source = "omics",
    top_n = 500,
    action = "add",
    family = "binomial"
)
```

    ## Log2-Transforming each assay in MultiAssayExperiment.

    ## Scaling each assay in MultiAssayExperiment.

    ## Running GLMs.

Now we can visualize these results with a manhattan plot.

``` r
# Plot the manhattan plot
plot_manhattan(
    exposomicset = expom,
    min_per_cat = 0,
    feature_col = "feature_clean",
    vars_to_label = c(
        "TC19001180.hg.1",
        "TC01000565.hg.1",
        "cg01937701",
        "hs_mibp_cadj_Log2"
    ),
    panel_sizes = c(1, 3, 1, 3, 1, 1, 1),
    facet_angle = 0
)
```

![Manhattan plot of omics-wide associations with asthma
status.](tidyexposomics_files/figure-html/plot-manhattan-1.png)

Manhattan plot of omics-wide associations with asthma status.

### Exposome Scores

------------------------------------------------------------------------

We can also calculate exposome scores, which are a summary measure of
exposure. The `run_exposome_score` function is used to calculate the
exposome score. The `exposure_cols` argument is used to set the columns
to use for the exposome score. The `score_type` argument is used to set
the type of score to calculate. Here we could use:

- `median`: Calculates the median of the exposure variables.

- `mean`: Calculates the mean of the exposure variables.

- `sum`: Calculates the sum of the exposure variables.

- `pca`: Calculates the first principal component of the exposure
  variables.

- `irt`: Uses Item Response Theory to calculate the exposome score.

- `quantile`: Calculates the quantile of the exposure variables.

- `var`: Calculates the variance of the exposure variables.

The `score_column_name` argument is used to set the name of the column
to store the exposome score in. Here we will define a score for aerosols
using a variety of different methods and demonstrate their use in
association with asthma status.

``` r
# determine which aerosol variables to use
aerosols <- c("h_pm25_ratio_preg_None", "h_pm10_ratio_preg_None")

# Create exposome scores
expom <- expom |>
    run_exposome_score(
        exposure_cols = aerosols,
        score_type = "median",
        score_column_name = "exposome_median_score"
    ) |>
    run_exposome_score(
        exposure_cols = aerosols,
        score_type = "pca",
        score_column_name = "exposome_pca_score"
    ) |>
    run_exposome_score(
        exposure_cols = aerosols,
        score_type = "irt",
        score_column_name = "exposome_irt_score"
    ) |>
    run_exposome_score(
        exposure_cols = aerosols,
        score_type = "quantile",
        score_column_name = "exposome_quantile_score"
    ) |>
    run_exposome_score(
        exposure_cols = aerosols,
        score_type = "var",
        score_column_name = "exposome_var_score"
    )
```

    ## Extracting exposure data...
    ## Extracting exposure data...
    ## Extracting exposure data...
    ## Extracting exposure data...
    ## Extracting exposure data...

    ## Calculating median exposure scores...

    ## Calculating PCA exposure scores...

    ## Calculating IRT exposure scores...

    ## Registered S3 method overwritten by 'vegan':
    ##   method     from      
    ##   rev.hclust dendextend

    ## Warning: EM cycles terminated after 500 iterations.

    ## Calculating quantile exposure scores...

    ## Calculating variance exposure scores...

We can then associate these exposome scores with asthma status using the
`run_association` function, just like we did before. However, this time
we specify our `feature_set` to be the exposome scores we just
calculated.

``` r
# Associate exposome scores with outcome
expom <- run_association(
    exposomicset = expom,
    outcome = "hs_asthma",
    source = "exposures",
    feature_set = c(
        "exposome_median_score",
        "exposome_pca_score",
        "exposome_irt_score",
        "exposome_quantile_score",
        "exposome_var_score"
    ),
    action = "add",
    family = "binomial"
)
```

    ## Running GLMs.

``` r
# Plot the association forest plot
plot_association(
    exposomicset = expom,
    source = "exposures",
    terms = c(
        "exposome_median_score",
        "exposome_pca_score",
        "exposome_irt_score",
        "exposome_quantile_score",
        "exposome_var_score"
    ),
    filter_col = "p.value",
    filter_thresh = 0.05,
    r2_col = "r2"
)
```

![Associations of aerosol exposome scores with asthma status. The
variance-based score has the strongest association with asthma
status.](tidyexposomics_files/figure-html/plot-exposome-scores-1.png)

Associations of aerosol exposome scores with asthma status. The
variance-based score has the strongest association with asthma status.

## Differential Abundance

------------------------------------------------------------------------

### Differential Abundance

We provide functionality to test for differentially abundant features
associated with an outcome across multiple omics layers. This is done
using the `run_differential_abundance` function, which fits a model
defined by the user (using the `formula` argument) and supports several
methods. Here we apply the `limma_trend` method, a widely used approach
for analyzing omics data. Users can also specify how features are scaled
(e.g. none, quantile, TMM) before fitting.

``` r
# Run differential abundance analysis
expom <- run_differential_abundance(
    exposomicset = expom,
    formula = ~hs_asthma,
    method = "limma_trend",
    scaling_method = "none",
    action = "add"
)
```

    ## Running differential abundance testing.

    ## Processing assay: Gene Expression

    ## Processing assay: Methylation

    ## Differential abundance testing completed.

We can summarize the results of the differential abundance analysis with
a volcano plot, which highlights features with a high log fold change
and that are statistically significant. The `plot_volcano` function
generates this visualization, with options to set thresholds for
p-values and log fold changes, and to label a subset of top-ranked
features. In this example, we use the `feature_clean` column to display
interpretable feature names.

**Note:** we set the `pval_col` to `P.Value` for the purposes of this
example, but we recommend keeping the default of `adj.P.Val` to use the
adjusted p-values.

``` r
# Plot the volcano plot
plot_volcano(
    exposomicset = expom,
    top_n_label = 3,
    feature_col = "feature_clean",
    logFC_thresh = log2(1),
    pval_thresh = 0.05,
    pval_col = "P.Value",
    logFC_col = "logFC",
    nrow = 1
)
```

![Volcano plot of differentially abundant features across omics
layers.](tidyexposomics_files/figure-html/volcano-plot-1.png)

Volcano plot of differentially abundant features across omics layers.

### Sensitivity Analysis

Depending on pre-processing steps, the results of the differential
abundance analysis may vary. The `sensitivity_analysis` function is used
to perform a sensitivity analysis to determine the robustness of the
results. Here we determine if a feature is still differentially abundant
if different minimum values, proportions, scaling methods are used, the
inclusion of covariates, and after bootstrapping. We then define a
stability score based on the number of times a feature is found to be
differentially abundant under different conditions as well as the
consistency of the effect size:

Where:

\\ Stability\\ Score = \frac{\sum_i{(p_i \< \alpha)}}{N} \* \frac{1}{1 +
\frac{\sigma\_{\beta}}{\mu\_{\|\beta\|}}}\\

Where:

- \\p_i\\ is the p-value for the \\i^{th}\\ test

- \\\alpha\\ is the significance threshold

- \\\beta\\ is the log fold change

- \\N\\ is the number of tests

- \\\sigma\_{\beta}\\ is the standard deviation of the effect size
  estimates

- \\\mu\_{\|\beta\|}\\ is the mean of the absolute value of the effect
  size estimates.

The first term captures the proportion of tests that are significant,
while the second term captures the consistency of the effect size
estimates. The stability score ranges from 0 to 1.

A stability score of 1 indicates that the feature is always found to be
differentially abundant, while a stability score of 0 indicates that the
feature is never found to be differentially abundant. Besides these, we
provide other score metrics as well:

- `presence_rate`: Proportion of runs in which the feature’s p-value is
  below the specified threshold (selection frequency).

- `effect_consistency`: Inverse of the coefficient of variation of log
  fold-changes; measures effect size stability across runs.

- `stability_score`: Hybrid score combining presence_rate and
  effect_consistency, capturing reproducibility and signal strength.

- `mean_log_p`: Average of negative log-transformed p-values; represents
  overall statistical signal strength.

- `logp_weighted_score`: Product of `mean_log_p` and effect_consistency;
  highlights consistently strong features.

- `sd_logFC`: Standard deviation of log fold-change estimates;
  quantifies variability of effect sizes.

- `iqr_logFC`: Interquartile range of log fold-changes; provides a
  robust measure of effect size spread.

- `cv_logFC`: Coefficient of variation of log fold-changes; reflects
  relative variability of effect size.

- `sign_flip_freq`: Proportion of runs where the sign of the effect size
  differs from the overall average direction.

- `sd_log_p`: Standard deviation of log-transformed p-values; indicates
  variability in statistical signal.

``` r
# Perform Sensitivity Analysis
expom <- run_sensitivity_analysis(
    exposomicset = expom,
    base_formula = ~hs_asthma,
    methods = c("limma_trend"),
    scaling_methods = c("none"),
    pval_col = "P.Value",
    logfc_col = "logFC",
    logFC_threshold = log2(1),
    pval_threshold = 0.05,
    stability_metric = "stability_score",
    bootstrap_n = 10,
    action = "add"
)
```

    ## Running bootstrap iteration 1 of 10

    ## Running bootstrap iteration 2 of 10

    ## Running bootstrap iteration 3 of 10

    ## Running bootstrap iteration 4 of 10

    ## Running bootstrap iteration 5 of 10

    ## Running bootstrap iteration 6 of 10

    ## Running bootstrap iteration 7 of 10

    ## Running bootstrap iteration 8 of 10

    ## Running bootstrap iteration 9 of 10

    ## Running bootstrap iteration 10 of 10

    ## Computing feature stability across sensitivity conditions.

    ## Feature stability analysis completed.

    ## Number Of Features Above Threshold Of 0.32:

    ## ----------------------------------------

    ## Gene Expression: 11/471

    ## Methylation: 5/277

Now we can plot the distribution of stability scores and visualize how
many features have scores greater than a threshold of 0.34.

``` r
# Plot distriubtion of stability scores
plot_sensitivity_summary(
    exposomicset = expom,
    stability_score_thresh = 0.34,
    stability_metric = "stability_score"
)
```

    ## Number of Features with stability_score > 0.34:

    ## Gene Expression: 9 / 471

    ## Methylation: 4 / 277

    ## Picking joint bandwidth of 0.024

![Sensitivity summary highlighting robust features by stability
score.](tidyexposomics_files/figure-html/plot-stability-score-1.png)

Sensitivity summary highlighting robust features by stability score.

## Multi-Omics Integration

------------------------------------------------------------------------

While differential abundance analysis per omics layer can deliver
insights into how each omic is associated with a particular outcome, we
may want to leverage methods which integrate multiple omics layers. The
`run_multiomics_integration` function is used to integrate multiple
omics layers. Here we use either `MCIA`, `RGCCA`, `MOFA`, or the
`DIABLO` method to integrate omics layers:

- `MCIA`: Multiple Co-inertia Analysis, a method that uses canonical
  correlation analysis to integrate multiple omics layers. We use the
  nipalsMCIA algorithm to compute the co-inertia scores from the
  [`nipalsMCIA`
  package](https://www.bioconductor.org/packages/release/bioc/html/nipalsMCIA.html)
  (nipalsMCIA).

- `RGCCA` : Generalized canonical correlation for flexible multi-block
  integration implemented using the [`RGCCA`
  package](https://rgcca-factory.github.io/RGCCA/) (Regularized and
  Sparse Generalized Ca…).

- `MOFA`: Multi-Omics Factor Analysis, a method that uses factor
  analysis to integrate multiple omics layers. MOFA is implemented using
  the [`MOFA2`
  package](https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html)
  (MOFA2).

- `DIABLO`: Supervised multi-block PLS for outcome-aligned latent
  factors implemented using the [`mixOmics`
  package](https://www.bioconductor.org/packages/devel/bioc/html/mixOmics.html)
  (mixOmics).

Here we are interested in integrating our omics layers with the end goal
of identifying multi-omics features that are associated with asthma
status. So, we will use the `DIABLO` method.

``` r
# Perform multi-omics integration
expom <- run_multiomics_integration(
    exposomicset = expom,
    method = "DIABLO",
    n_factors = 5,
    outcome = "hs_asthma",
    action = "add"
)
```

    ## Scaling each assay in MultiAssayExperiment.

    ## Running multi-omics integration using DIABLO...

    ## Applying DIABLO supervised integration.

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

We can then use `plot_factor_summary` to visualize which omics
contribute most to which factors.

``` r
# Plot factor summary
plot_factor_summary(
    exposomicset = expom,
    midpoint = 2.5
)
```

![Contribution of each omics layer to DIABLO latent
factors.](tidyexposomics_files/figure-html/factor-summary-1.png)

Contribution of each omics layer to DIABLO latent factors.

Here we see that these factors are largely driven by features in the
gene expression data.

### Factor Analysis

These methods are designed to identify factors that we can then
associate with an outcome variable. Here we will use the
`run_association` function to identify factors that are associated with
asthma status.

``` r
# Identify factors that are associated with the outcome
expom <- run_association(
    exposomicset = expom,
    source = "factors",
    outcome = "hs_asthma",
    feature_set = exp_vars,
    action = "add",
    family = "binomial"
)
```

    ## Running GLMs.

Now let’s see if any of our factors are associated with asthma status.

``` r
# Plot the forest plot
plot_association(
    exposomicset = expom,
    source = "factors",
    filter_col = "p_adjust",
    filter_thresh = 0.05,
    r2_col = "r2"
)
```

![Associations of latent factors with asthma after covariate
adjustment.](tidyexposomics_files/figure-html/plot-associate-factors-1.png)

Associations of latent factors with asthma after covariate adjustment.

We see that several factors are associated with our outcome. Factors
have loading scores which indicate the strength of the association
between the factor and the features. Here we can extract the top
features, those in the 70th percentile, per omic, associated with our
factors of interest. We set the `pval_col` and `pval_thresh` to filter
our association results to grab the factors that pass our thresholds.
However, we could specify specific factors with the `factors` argument
too.

``` r
# Extract top features that contribute to a factor
expom <- extract_top_factor_features(
    exposomicset = expom,
    method = "percentile",
    pval_col = "p_adjust",
    pval_thresh = 0.05,
    percentile = 0.7,
    action = "add"
)
```

    ## Extracting top contributing features for specified factors.

    ## Using DIABLO loadings.

    ## Applying percentile-based filtering (>70%).

    ## Selected 589 features contributing to specified factors.

We can visualize the top features associated with each factor using the
`plot_top_factor_features` function. The `top_n` argument is used to set
the number of top features to plot. The factors argument is used to set
the factors to plot.

``` r
# Plot the top factor features
plot_top_factor_features(
    exposomicset = expom,
    top_n = 5,
    feature_col = "feature_clean"
)
```

![Top-loading features per factor and omics
layer.](tidyexposomics_files/figure-html/plot-top-features-1.png)

Top-loading features per factor and omics layer.

## Exposure-Omics Association

------------------------------------------------------------------------

### Exposure-Omics Association

Now we have the option to correlate either the top factor features,
differentially abundant features, or user-specified omics features (by
using a variable map, a data frame with two columns, `exp_name` for the
name of the omics assay, and `variable` for the name of the molecular
feature) with exposures.

Here we will correlate features driving multiple latent factors with
exposures. To grab these features, we will grab them from the metadata
in our `MultiAssayExperiment` object. The `correlation_cutoff` is used
to set the minimum correlation score, while the `pval_cutoff` is used to
set the maximum p-value for the association.

``` r
# Grab top common factor features and ensure
# feature is renamed to variable for the variable_map
top_factor_features <- expom |>
    extract_results(result = "multiomics_integration") |>
    pluck("top_factor_features") |>
    dplyr::select(variable = feature, exp_name)

# Correlate top factor features with exposures
# Perform correlation analysis between
# factor features and exposures
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "omics",
    variable_map = top_factor_features,
    exposure_cols = exp_vars,
    action = "add",
    correlation_cutoff = 0.3,
    pval_cutoff = 0.05,
    cor_pval_column = "p.value"
)
```

    ## Log2-Transforming each assay in MultiAssayExperiment.

``` r
# Perform correlation analysis between factor features
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "omics",
    variable_map = top_factor_features,
    exposure_cols = exp_vars,
    feature_cors = TRUE,
    action = "add",
    correlation_cutoff = 0.3,
    pval_cutoff = 0.05,
    cor_pval_column = "p.value"
)
```

    ## Log2-Transforming each assay in MultiAssayExperiment.

We can plot the results of the exposure-omics association analysis using
the `plot_correlation_summary` function. Here we set the mode to
`summary`, which will plot the number of associations per exposure and
feature type. The `feature_type` argument is used to set the type of
features to plot. Here we use `omics`, which are the top factor
features.

``` r
# Plot correlation summary
plot_correlation_summary(
    exposomicset = expom,
    mode = "summary",
    feature_type = "omics"
)
```

![Summary of exposure-omics and feature-feature
associations.](tidyexposomics_files/figure-html/exp-omic-corr-summary-1.png)

Summary of exposure-omics and feature-feature associations.

Here we note that gene expression features have a greater number of
associations with exposures. Among exposures, we see that aerosols have
the greatest number of associations with features, which may reflect the
larger number of aerosol variables.

It may also be useful to identify which exposures are correlated with
similar molecular features. We can do this with the
`plot_circos_correlation` function. This function will plot a circos
plot of the exposures and their shared features. The `feature_type`
argument is used to set the correlation results to use for the analysis.
Here we use the `omics` feature set, which plots the feature-exposure
correlation information. The cutoff argument is used to set the minimum
number of shared features to plot. Here we set the `shared_cutoff` to 1,
which means that only exposures that share at least one feature will be
plotted.

``` r
# Plot the circos plot of shared correlated features
plot_circos_correlation(
    exposomicset = expom,
    feature_type = "omics",
    shared_cutoff = 1,
    midpoint = 1.5
)
```

![Circos of exposures sharing correlated molecular
features.](tidyexposomics_files/figure-html/exposure-shared-corr-1.png)

Circos of exposures sharing correlated molecular features.

Here we see that particulate matter exposures are correlated with many
of the same features.

### Network Analysis

The `run_create_network` function is used to create a network of
exposures and omics features. The correlation results are used to create
the network. The `feature_type` argument is used to set the type of
features to use for the network. We could choose any of the following:

- `degs`: Differentially abundant features correlated with exposures.

- `factors`: Factor features correlated with exposures.

- `omics`: User specified omics features correlated with exposures.

Here we create networks based on two correlation tables generated above:

- `omics`: Correlations between omics and exposure variables.

- `omics_feature_cor`: Correlation just between the omics features.

This allows us to quantify a bipartite graph between exposures and omics
features and a graph just between the omics features.

``` r
# Create omics within feature correlation network
expom <- run_create_network(
    exposomicset = expom,
    feature_type = "omics_feature_cor",
    action = "add"
)
```

    ## Creating network from correlation results.

    ## Network added to metadata as: network_omics_feature_cor

``` r
# Create omics-exposure correlation network
expom <- run_create_network(
    exposomicset = expom,
    feature_type = "omics",
    action = "add"
)
```

    ## Creating network from correlation results.

    ## Network added to metadata as: network_omics

To plot the network we can use the `plot_network` function. The network
argument is used to set the type of network to plot. The `top_n_nodes`
argument is used to set the number of nodes to plot. The
`node_color_var` argument is used to set the variable to color the nodes
by. The `label` argument is used to set whether or not to label the
nodes. We can label the top n nodes, where the default is 5 nodes based
on centrality. Additionally, we can choose to label certain nodes using
the `nodes_to_label` argument. We can also include the network
statistics using the `include_stats` argument.

``` r
# Plot the exposure-omics network
# Setting the seed so that the graph layout is consistent 
# given that layouts can have random starting points
set.seed(42)
plot_network(
    exposomicset = expom,
    network = "omics",
    top_n_nodes = 50,
    include_stats = TRUE,
    cor_thresh = 0.2,
    node_color_var = "group",
    label = TRUE,
    label_top_n = 5
)
```

    ## Extracting graph.

![Bipartite network of exposures and correlated omics features where
node labels reflect most central
nodes.](tidyexposomics_files/figure-html/plot-network-1.png)

Bipartite network of exposures and correlated omics features where node
labels reflect most central nodes.

Here we note that the network is highly connected, with particulate
matter and methyl paraben exposures having the most connections to omics
features. To briefly explain the metrics mentioned:

- **Nodes**: Number of nodes included in the graph.

- **Edges**: Number of connections between nodes in the graph.

- **Components**: Number of disconnected subgraphs.

- **Diameter**: The longest shortest path between any two nodes.

- **Mean Distance**: The average shortest path between any two nodes.

- **Modularity**: A measure of the division in the graph, where higher
  modularity indicates that nodes within the same community are more
  densely connected compared to nodes in different communities.

### Exposure-Omics Impact Network Analysis

The bipartite graph above describes how exposures are associated with
omics features. However, this does not describe if exposures are
associated with features that are more central in the feature-only graph
(i.e. are exposures associated with features that are more central or
peripheral to the graph?). To answer this question, we can use the
`run_exposure_impact` function to calculate centrality metrics for the
omics features associated with a given exposure. Centrality is computed
on the feature-only graph. Centrality metrics are used to identify
features that are more central to the network and may be more important
in the context of the exposure-omics relationships. The centrality
metrics included are:

- **Mean Betweenness Centrality**: Measures how often a node lies on the
  shortest path between other nodes where high values indicate potential
  intermediates.

- **Mean Closeness Centrality**: Measures how close a node is to all
  others in the network where high values indicate faster access to all
  nodes.

- **Mean Degree**: The average number of direct connections or edges a
  node has.

- **Mean Eigenvector Centrality**: Measures how well-connected a node is
  and how well-connected its neighbors are, where high values suggest
  influence in a well-connected cluster.

``` r
# Run exposure-omics impact analysis
expom <- run_exposure_impact(
    exposomicset = expom,
    feature_type = "omics"
)
```

To plot the results of the exposure-omics impact analysis, we can use
the `plot_exposure_impact` function. We set `feature_type` to `omics`,
which means that we are plotting the network metrics of the omics
features correlated with exposures. The `min_per_group` argument is used
to set the minimum number of features per group to plot. This is useful
for filtering out exposures that do not have enough associated features.

``` r
# Plot the exposure impact summary
plot_exposure_impact(
    exposomicset = expom,
    feature_type = "omics",
    min_per_group = 10
)
```

![Network centrality metrics of omics features associated with each
exposure.](tidyexposomics_files/figure-html/plot-exposure-impact-1.png)

Network centrality metrics of omics features associated with each
exposure.

Here we see that child PFNA (`hs_pfna_c_Log2`) exposure is associated
with omics features that are central to our network, despite not having
the greatest number of features associated with it. This suggests that
childhood PFNA exposure may have a significant impact on the omics
features, even if it is not associated with as many features as some of
the other exposures.

## Enrichment Analysis

------------------------------------------------------------------------

Enrichment analysis tests whether a set of molecular features
(e.g. differentially abundant genes, metabolites, etc.) is
over-represented in a predefined biological process. The benefit of
grouping our exposures into categories is that we can now determine how
broad categories of exposures are tied to biological processes. The
`run_enrichment` function can perform enrichment analysis on the
following feature types:

- `degs`: Differentially abundant features.

- `degs_robust`: Robust differentially abundant features from the
  sensitivity analysis.

- `omics`: User chosen features.

- `factor_features`: Multi-omics factor features either from
  `factor_type = “common_top_factor_features”` or
  `“top_factor_features”`.

- `degs_cor`: Differentially abundant features correlated with a set of
  exposures.

- `omics_cor`: User chosen features correlated with a set of exposures.

- `factor_features_cor`: Multi-omics factor features correlated with a
  set of exposures.

Here we will run enrichment analysis on factor features correlated with
exposures. We specify `feature_col` to represent the column in our
feature metadata with IDs that can be mapped (i.e. gene names). We will
be performing Gene ontology enrichment powered by the [`fenr`
package](https://www.bioconductor.org/packages/release/bioc/html/fenr.html)
(Fenr, 2025). Note that we specify a `clustering_approach.` This will
cluster our enrichment terms by the molecular feature overlap.

``` r
# Run enrichment analysis on factor features correlated with exposures
expom <- run_enrichment(
    exposomicset = expom,
    feature_type = c("omics_cor"),
    feature_col = "feature_clean",
    db = c("GO"),
    species = "goa_human",
    fenr_col = "gene_symbol",
    padj_method = "none",
    pval_thresh = 0.1,
    min_set = 1,
    max_set = 800,
    clustering_approach = "diana",
    action = "add"
)
```

    ## Determining Number of GO Term Clusters...

    ## Optimal number of clusters for samples: 6

### Enrichment Visualizations

To visualize our enrichment results we provide several options:

- dot\`plot: A dot plot showing the top enriched terms. The size of the
  dots represents the number of features associated with the term, while
  the color represents the significance of the term.

- `cnet`: A network plot showing the relationship between features and
  enriched terms.

- `network`: A network plot showing the relationship between enriched
  terms.

- `heatmap`: A heatmap showing the relationship between features and
  enriched terms.

- `summary`: A summary figure of the enrichment results.

#### Enrichment Summary

To summarize the enrichment results, we can use the `plot_enrichment`
function with the `plot_type` argument set to `summary.` This will plot
a summary of the enrichment results, showing:

- The number of exposure categories per enrichment term group.

- The number of features driving the enrichment term group.

- A p-value distribution of the enrichment term group.

- The number of terms in the enrichment term group.

- The total number of terms per experiment name.

- The overlap in enrichment terms between experiments (i.e. between gene
  expression and methylation).

``` r
# Plot the summary diagram
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    plot_type = "summary"
)
```

    ## Picking joint bandwidth of 0.303

![Summary of enriched GO terms grouped by overlap and exposure
category.](tidyexposomics_files/figure-html/enrichment-summary-1.png)

Summary of enriched GO terms grouped by overlap and exposure category.

Here we see that it is just the features associated with “polyatomic
entity” exposures that seem to be enriched. Additionally, there appears
to be no overlap in terms between methylation and gene expression
results.

#### DotPlot

By setting the `plot_type` to dotplot we can create a dotplot to show
which omics and which exposure categories are associated with which
terms. By specifying the `top_n_genes` we can add the most frequent
features in that particular enrichment term group.

``` r
# Plot a dotplot of terms
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    plot_type = "dotplot",
    top_n = 15,
    add_top_genes = TRUE,
    top_n_genes = 5
)
```

![Dotplot of top enriched GO terms by omics layer and exposure
category.](tidyexposomics_files/figure-html/enrichment-dotplot-1.png)

Dotplot of top enriched GO terms by omics layer and exposure category.

#### Term Network Plot

We can set the `plot_type` to `network` to understand how our enrichment
terms are individually connected.

``` r
# Plot the term network plot
# Setting a seed so that the plot layout is consistent
set.seed(42)
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    plot_type = "network",
    label_top_n = 1
)
```

![Network of enriched GO terms connected by shared
genes.](tidyexposomics_files/figure-html/enrichment-term-network-1.png)

Network of enriched GO terms connected by shared genes.

At the individual term level, we see that they differ by omics layer,
with the gene expression driving terms related to protein binding and
the methylation data driving terms related to DNA binding and
transcription.

#### Heatmap

Setting the `plot_type` to `heatmap` can help us understand which genes
are driving the enrichment terms. We have the additional benefit of
being able to color our tiles by the Log_2_Fold Change from our
differential abundance testing. Here we will examine group 1, given it
seems to be driven by the most terms and multiple omics layers.

``` r
# Plot a heatmap of genes and corresponding GO terms
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    go_groups = "Group 1",
    plot_type = "heatmap",
    heatmap_fill = TRUE,
    feature_col = "feature_clean"
)
```

![Heatmap of genes driving enriched GO terms (Group 1) with log2
fold-change
overlay.](tidyexposomics_files/figure-html/enrichment-heatmap-1.png)

Heatmap of genes driving enriched GO terms (Group 1) with log2
fold-change overlay.

Here we see that methylation features have the strongest fold change
with respect to asthma status. Interestingly, we see that gene
expression of *GZMH*, a potent immune effector molecule, is upregulated
in asthma and correlated with exposures in the polyatomic entity
category.

#### Cnet Plot

Another way to visualize this information is with the `cnet` plot, where
the enrichment terms are connected to the genes driving them.

``` r
# Plot the gene-term network
# Setting a seed so that the plot layout is consistent
set.seed(42)
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    go_groups = "Group 1",
    plot_type = "cnet",
    feature_col = "feature_clean"
)
```

![Cnet plot linking enriched terms to contributing genes (Group
1).](tidyexposomics_files/figure-html/enrichment-cnetplot-1.png)

Cnet plot linking enriched terms to contributing genes (Group 1).

## Custom Analysis

We provide functionality to access the underlying data in the
MultiAssayExperiment object and to construct tibbles for your own
analysis:

- `pivot_sample`: Pivot the sample data to a tibble with samples as rows
  and exposures as columns.

- `pivot_feature`: Pivot the feature metadata to a tibble with features
  as rows and feature metadata as columns.

- `pivot_exp`: Pivot the sample and experiment assay data to a tibble
  with samples as rows and sample metadata as columns. Additionally,
  there will be a column for values for specified features in specified
  assays.

### Pivot Sample

Let’s check out the `pivot_sample` function. This function pivots the
sample data to a tibble with samples as rows and exposures as columns.

``` r
# Pivot sample data to a tibble
expom |>
    pivot_sample() |>
    head()
```

    ## # A tibble: 6 × 73
    ##   .sample e3_sex_None hs_child_age_None h_age_None h_mbmi_None e3_gac_None
    ##   <chr>   <fct>                   <dbl>      <dbl>       <dbl>       <dbl>
    ## 1 s812    female                   6.84       35          33.2        41.4
    ## 2 s296    female                   6.75       32.8        25.1        39.1
    ## 3 s242    female                   6.97       33          24.6        40.4
    ## 4 s376    male                     6.35       27.1        28.3        37  
    ## 5 s1183   female                   6.74       32.7        33.3        38.1
    ## 6 s98     female                   6.93       23          31.2        39  
    ## # ℹ 67 more variables: h_native_None <fct>, h_parity_None <fct>,
    ## #   hs_asthma <dbl>, hs_zbmi_who <dbl>, h_pm10_ratio_preg_None <dbl>,
    ## #   h_pm25_ratio_preg_None <dbl>, hs_pm10_yr_hs_h_None <dbl>,
    ## #   hs_pm25_yr_hs_h_None <dbl>, hs_pb_c_Log2 <dbl>, hs_pfhxs_m_Log2 <dbl>,
    ## #   hs_pfna_c_Log2 <dbl>, hs_bpa_madj_Log2 <dbl>, hs_mibp_cadj_Log2 <dbl>,
    ## #   PC1 <dbl>, PC2 <dbl>, PC3 <dbl>, PC4 <dbl>, PC5 <dbl>, PC6 <dbl>,
    ## #   PC7 <dbl>, PC8 <dbl>, PC9 <dbl>, PC10 <dbl>, PC11 <dbl>, PC12 <dbl>, …

We could use this functionality to count the number of asthmatics per
sex:

``` r
# Count the number of asthma cases by sex status
expom |>
    pivot_sample() |>
    group_by(hs_asthma, e3_sex_None) |>
    reframe(n = n())
```

    ## # A tibble: 4 × 3
    ##   hs_asthma e3_sex_None     n
    ##       <dbl> <fct>       <int>
    ## 1         0 female         23
    ## 2         0 male           15
    ## 3         1 female          5
    ## 4         1 male            4

### Pivot Feature

The `pivot_feature` function pivots the feature metadata to a tibble
with features as rows and feature metadata as columns. This can be
useful for exploring the feature metadata in a more flexible way.

``` r
# Pivot feature data to a tibble
expom |>
    pivot_feature() |>
    head()
```

    ## # A tibble: 6 × 42
    ##   .exp_name     .feature transcript_cluster_id probeset_id seqname strand  start
    ##   <chr>         <chr>    <chr>                 <chr>       <chr>   <chr>   <int>
    ## 1 Gene Express… TC02004… TC02004973.hg.1       TC02004973… chr2    +      8.99e7
    ## 2 Gene Express… TC02004… TC02004949.hg.1       TC02004949… chr2    -      8.96e7
    ## 3 Gene Express… TC14002… TC14002235.hg.1       TC14002235… chr14   -      1.06e8
    ## 4 Gene Express… TC15000… TC15000863.hg.1       TC15000863… chr15   +      9.09e7
    ## 5 Gene Express… TC19001… TC19001585.hg.1       TC19001585… chr19   -      4.39e7
    ## 6 Gene Express… TC17000… TC17000226.hg.1       TC17000226… chr17   +      1.90e7
    ## # ℹ 35 more variables: stop <int>, total_probes <int>, gene_assignment <chr>,
    ## #   mrna_assignment <chr>, notes <chr>, phase <chr>, TC_size <dbl>,
    ## #   TSS_Affy <int>, EntrezeGeneID_Affy <chr>, GeneSymbol_Affy <chr>,
    ## #   GeneSymbolDB <chr>, GeneSymbolDB2 <chr>, mrna_ID <chr>, mrna_DB <chr>,
    ## #   mrna_N <int>, notesYN <chr>, geneYN <chr>, genes_N <dbl>, CallRate <dbl>,
    ## #   fil1 <chr>, feature_clean <chr>, Forward_Sequence <chr>, SourceSeq <chr>,
    ## #   Random_Loci <chr>, Methyl27_Loci <chr>, UCSC_RefGene_Name <chr>, …

We could use this functionality to count the number of features per
omics layer, or to filter features based on their metadata. For example,
we can count the number of features per omics layer:

``` r
# Count the number of features per omic layer
expom |>
    pivot_feature() |>
    group_by(.exp_name) |>
    summarise(n = n())
```

    ## # A tibble: 2 × 2
    ##   .exp_name           n
    ##   <chr>           <int>
    ## 1 Gene Expression   471
    ## 2 Methylation       277

### Pivot Experiment

Now if we want to grab assay data from a particular experiment, we can
do that with the `pivot_exp.` Let’s try grabbing assay values for the
`TC01004453.hg.1` probe (probe for *IL23R*) in the `Gene Expression`
experiment:

``` r
# Pivot experiment data to a tibble
expom |>
    pivot_exp(
        exp_name = "Gene Expression",
        features = "TC01004453.hg.1"
    ) |>
    head()
```

    ## # A tibble: 6 × 102
    ##   exp_name     .feature .sample counts transcript_cluster_id probeset_id seqname
    ##   <chr>        <chr>    <chr>    <dbl> <chr>                 <chr>       <chr>  
    ## 1 Gene Expres… TC01004… s812     1.09  TC01004453.hg.1       TC01004453… chr1   
    ## 2 Gene Expres… TC01004… s296     1.71  TC01004453.hg.1       TC01004453… chr1   
    ## 3 Gene Expres… TC01004… s242     1.19  TC01004453.hg.1       TC01004453… chr1   
    ## 4 Gene Expres… TC01004… s376     1.21  TC01004453.hg.1       TC01004453… chr1   
    ## 5 Gene Expres… TC01004… s1183    0.803 TC01004453.hg.1       TC01004453… chr1   
    ## 6 Gene Expres… TC01004… s98      0.835 TC01004453.hg.1       TC01004453… chr1   
    ## # ℹ 95 more variables: strand <chr>, start <int>, stop <int>,
    ## #   total_probes <int>, gene_assignment <chr>, mrna_assignment <chr>,
    ## #   notes <chr>, phase <chr>, TC_size <dbl>, TSS_Affy <int>,
    ## #   EntrezeGeneID_Affy <chr>, GeneSymbol_Affy <chr>, GeneSymbolDB <chr>,
    ## #   GeneSymbolDB2 <chr>, mrna_ID <chr>, mrna_DB <chr>, mrna_N <int>,
    ## #   notesYN <chr>, geneYN <chr>, genes_N <dbl>, CallRate <dbl>, fil1 <chr>,
    ## #   feature_clean <chr>, e3_sex_None <fct>, hs_child_age_None <dbl>, …

We can use this functionality to create custom plots or analyses based
on the exposure and feature data. For example, we can plot the
expression levels of *IL23R* by asthma status:

``` r
# Plot expression of IL23R
expom |>
    pivot_exp(
        exp_name = "Gene Expression",
        features = "TC01004453.hg.1"
    ) |>
    ggplot(aes(
        x = as.character(hs_asthma),
        y = log2(counts + 1),
        color = as.character(hs_asthma),
        fill = as.character(hs_asthma)
    )) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(alpha = 0.1) +
    ggpubr::geom_pwc(label = "{p.adj.format}{p.adj.signif}") +
    theme_minimal() +
    ggpubr::rotate_x_text(angle = 45) +
    ggsci::scale_color_cosmic() +
    ggsci::scale_fill_cosmic() +
    labs(
        x = "",
        y = expression(Log[2] * "Abd."),
        fill = "Asthma Status",
        color = "Asthma Status"
    )
```

![IL23R gene expression levels (probe TC01004453.hg.1) by asthma
status.](tidyexposomics_files/figure-html/pivot-exp-example-1.png)

IL23R gene expression levels (probe TC01004453.hg.1) by asthma status.

## Pipeline Summary

To summarize the steps we have taken in this analysis, we can use the
`run_pipeline_summary` function. This function will provide a summary of
the steps taken in the analysis. We can set `console_print` to `TRUE` to
print the summary to the console. Setting `include_notes` to `TRUE` will
include notes on the steps taken in the analysis.

``` r
# Run the pipeline summary
expom |>
    run_pipeline_summary(console_print = TRUE, include_notes = TRUE)
```

    ## 1. run_impute_missing - 
    ## 2. filter_omics_Methylation - Filtered omics features from 'Methylation'
    ##         Using method = 'variance': 223 removed of 500 (44.6%).
    ## 3. filter_omics_Gene Expression - Filtered omics features from 'Gene Expression'
    ##         Using method = 'expression': 29 removed of 500 (5.8%).
    ## 4. run_normality_check - Assessed normality of 15 numeric exposure variables. 9 were normally distributed (p > 0.05), 6 were not.
    ## 5. transform_exposure - Applied 'boxcox_best' transformation to 9 exposure variables. 5 passed normality (Shapiro-Wilk p > 0.05, 55.6%).
    ## 6. run_pca - Outliers:  s1231
    ## 7. filter_sample_outliers - Outliers:  s1231
    ## 8. run_correlation_pcs - Correlated pcs features with exposures.
    ## 9. run_cluster_samples - Optimal number of clusters for samples: 2
    ## 10. run_correlation_exposures - Correlated exposures features with exposures.
    ## 11. run_association - Performed association analysis using source: exposures
    ## 12. run_association - Performed association analysis using source: exposures
    ## 13. run_exposome_score_exposome_median_score - Exposome score computed using method: 'median'
    ## 14. run_exposome_score_exposome_pca_score - Exposome score computed using method: 'pca'
    ## 15. run_exposome_score_exposome_irt_score - Exposome score computed using method: 'irt'
    ## 16. run_exposome_score_exposome_quantile_score - Exposome score computed using method: 'quantile'
    ## 17. run_exposome_score_exposome_var_score - Exposome score computed using method: 'var'
    ## 18. run_association - Performed association analysis using source: exposures
    ## 19. run_differential_abundance - Performed differential abundance analysis across all assays.
    ## 20. run_sensitivity_analysis - Ran sensitivity analysis across 10 bootstrap iterations and 1 methods and 1 scaling strategies. No covariates were removed from the base model.
    ## 21. run_multiomics_integration - Performed multi-omics integration using DIABLO with 5 latent factors. Scaling was enabled.
    ## 22. run_association - Performed association analysis using source: exposures
    ## 23. extract_top_factor_features - Selected 589 features contributing to specified factors.
    ## 24. run_correlation_omics - Correlated omics features with exposures.
    ## 25. run_correlation_omics - Correlated omics features with exposures.
    ## 26. run_create_network - Created undirected network from correlation results for 'omics_feature_cor'.
    ## 27. run_create_network - Created undirected network from correlation results for 'omics_feature_cor'.
    ## 28. run_exposure_impact - Computed exposure impact using omics correlation network.
    ## 29. run_enrichment - Performed GO enrichment on omics_cor features.

### Saving Results

We can export the results in our `MultiAssayExperiment` to an Excel
spreadsheet using the `extract_results_excel` function. Here we add all
of our results to the Excel file, but we can choose certain results by
changing the `result_types` argument.

``` r
# Save results
extract_results_excel(
    exposomicset = expom,
    file = tempfile(),
    result_types = "all"
)
```

    ## Writing Correlation Results.

    ## Writing Association Results.

    ## Writing Mixture Analysis Results.

    ## Writing Differential Abundance Results.

    ## Writing Sensitivity Analysis Results.

    ## Writing Multiomics Integration Results.

    ## Writing Network Impact Results.

    ## Writing Enrichment Results.

    ## Writing Pipeline Summary.

    ## Writing Exposure Summary Results.

    ## Results written to: C:\Users\Jason\AppData\Local\Temp\RtmpawvHtg\file5d90418c5f06

### Session Info

``` r
sessionInfo()
```

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] tidyexposomics_0.99.14      MultiAssayExperiment_1.36.1
    ##  [3] SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ##  [5] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    ##  [7] IRanges_2.44.0              S4Vectors_0.48.0           
    ##  [9] BiocGenerics_0.56.0         generics_0.1.4             
    ## [11] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ## [13] lubridate_1.9.5             forcats_1.0.1              
    ## [15] stringr_1.6.0               dplyr_1.1.4                
    ## [17] purrr_1.2.0                 readr_2.2.0                
    ## [19] tidyr_1.3.2                 tibble_3.3.0               
    ## [21] ggplot2_4.0.2               tidyverse_2.0.0            
    ## [23] BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] R.methodsS3_1.8.2     vroom_1.7.0           rARPACK_0.11-0       
    ##   [4] tidybulk_2.0.1        nnet_7.3-20           DT_0.34.0            
    ##   [7] vctrs_0.6.5           digest_0.6.39         png_0.1-8            
    ##  [10] corpcor_1.6.10        shape_1.4.6.1         proxy_0.4-29         
    ##  [13] BiocBaseUtils_1.12.0  ggrepel_0.9.7         parallelly_1.46.1    
    ##  [16] permute_0.9-10        magick_2.9.1          MASS_7.3-65          
    ##  [19] pkgdown_2.2.0         reshape2_1.4.5        httpuv_1.6.16        
    ##  [22] foreach_1.5.2         withr_3.0.2           xfun_0.54            
    ##  [25] ggpubr_0.6.3          survival_3.8-3        doRNG_1.8.6.3        
    ##  [28] memoise_2.0.1         fenr_1.8.1            janitor_2.2.1        
    ##  [31] SimDesign_2.24        tidyHeatmap_1.13.1    ggsci_4.2.0          
    ##  [34] systemfonts_1.3.2     ragg_1.5.0            GlobalOptions_0.1.3  
    ##  [37] pbapply_1.7-4         R.oo_1.27.1           Formula_1.2-5        
    ##  [40] ellipse_0.5.0         promises_1.5.0        otel_0.2.0           
    ##  [43] httr_1.4.8            beepr_2.0             rstatix_0.7.3        
    ##  [46] globals_0.19.0        rstudioapi_0.18.0     missForest_1.6.1     
    ##  [49] base64enc_0.1-6       curl_7.0.0            ggraph_2.2.2         
    ##  [52] polyclip_1.10-7       randomForest_4.7-1.2  SparseArray_1.10.8   
    ##  [55] xtable_1.8-8          desc_1.4.3            doParallel_1.0.17    
    ##  [58] evaluate_1.0.5        S4Arrays_1.10.1       BiocFileCache_3.0.0  
    ##  [61] hms_1.1.4             bookdown_0.46         colorspace_2.1-2     
    ##  [64] filelock_1.0.3        magrittr_2.0.4        snakecase_0.11.1     
    ##  [67] later_1.4.8           viridis_0.6.5         lattice_0.22-7       
    ##  [70] future.apply_1.20.2   class_7.3-23          Hmisc_5.2-5          
    ##  [73] pillar_1.11.1         nlme_3.1-168          iterators_1.0.14     
    ##  [76] compiler_4.5.1        RSpectra_0.16-2       stringi_1.8.7        
    ##  [79] gower_1.0.2           dendextend_1.19.1     plyr_1.8.9           
    ##  [82] crayon_1.5.3          abind_1.4-8           naniar_1.1.0         
    ##  [85] mixOmics_6.34.0       locfit_1.5-9.12       graphlayouts_1.2.3   
    ##  [88] bit_4.6.0             codetools_0.2-20      textshaping_1.0.4    
    ##  [91] recipes_1.3.1         bslib_0.10.0          e1071_1.7-17         
    ##  [94] GetoptLong_1.1.0      mime_0.13             splines_4.5.1        
    ##  [97] circlize_0.4.17       Rcpp_1.1.1            dbplyr_2.5.2         
    ## [100] knitr_1.51            blob_1.3.0            utf8_1.2.6           
    ## [103] clue_0.3-67           itertools_0.1-3       fs_1.6.6             
    ## [106] listenv_0.10.0        checkmate_2.3.4       Rdpack_2.6.6         
    ## [109] densityClust_0.3.3    openxlsx_4.2.8.1      ggsignif_0.6.4       
    ## [112] Matrix_1.7-4          statmod_1.5.1         tzdb_0.5.0           
    ## [115] visdat_0.6.0          tweenr_2.0.3          pkgconfig_2.0.3      
    ## [118] tools_4.5.1           cachem_1.1.0          rbibutils_2.4.1      
    ## [121] RSQLite_2.4.6         viridisLite_0.4.3     rvest_1.0.5          
    ## [124] DBI_1.3.0             fastmap_1.2.0         rmarkdown_2.30       
    ## [127] scales_1.4.0          grid_4.5.1            audio_0.1-12         
    ## [130] broom_1.0.12          sass_0.4.10           patchwork_1.3.2      
    ## [133] FNN_1.1.4.1           BiocManager_1.30.27   carData_3.0-6        
    ## [136] rpart_4.1.24          farver_2.1.2          tidygraph_1.3.1      
    ## [139] mgcv_1.9-3            yaml_2.3.12           foreign_0.8-90       
    ## [142] cli_3.6.5             lifecycle_1.0.5       caret_7.0-1          
    ## [145] lava_1.8.2            sessioninfo_1.2.3     backports_1.5.0      
    ## [148] BiocParallel_1.44.0   timechange_0.4.0      gtable_0.3.6         
    ## [151] rjson_0.2.23          ggridges_0.5.7        progressr_0.18.0     
    ## [154] parallel_4.5.1        pROC_1.19.0.1         testthat_3.3.2       
    ## [157] limma_3.66.0          jsonlite_2.0.0        edgeR_4.8.2          
    ## [160] bit64_4.6.0-1         assertthat_0.2.1      brio_1.1.5           
    ## [163] Rtsne_0.17            vegan_2.7-3           zip_2.3.3            
    ## [166] ranger_0.18.0         GPArotation_2025.3-1  jquerylib_0.1.4      
    ## [169] R.utils_2.13.0        timeDate_4052.112     dcurver_0.9.3        
    ## [172] shiny_1.13.0          dynamicTreeCut_1.63-1 htmltools_0.5.9      
    ## [175] rappdirs_0.3.4        RGCCA_3.0.3           glue_1.8.0           
    ## [178] factoextra_2.0.0      ggvenn_0.1.19         httr2_1.2.2          
    ## [181] XVector_0.50.0        mirt_1.45.1           gridExtra_2.3        
    ## [184] igraph_2.2.2          R6_2.6.1              Deriv_4.2.0          
    ## [187] labeling_0.4.3        ggh4x_0.3.1           cluster_2.1.8.2      
    ## [190] rngtools_1.5.2        clipr_0.8.0           ipred_0.9-15         
    ## [193] DelayedArray_0.36.0   tidyselect_1.2.1      htmlTable_2.4.3      
    ## [196] ggforce_0.5.0         xml2_1.5.2            car_3.1-5            
    ## [199] future_1.69.0         ModelMetrics_1.2.2.2  S7_0.2.1             
    ## [202] data.table_1.18.2.1   htmlwidgets_1.6.4     ComplexHeatmap_2.26.1
    ## [205] RColorBrewer_1.1-3    rlang_1.1.7           ggnewscale_0.5.2     
    ## [208] Cairo_1.7-0           hardhat_1.4.2         prodlim_2025.04.28

### References

Chung, M. K., House, J. S., Akhtari, F. S., Makris, K. C., Langston, M.
A., Islam, K. T., Holmes, P., Chadeau-Hyam, M., Smirnov, A. I., Du, X.,
Thessen, A. E., Cui, Y., Zhang, K., Manrai, A. K., Motsinger-Reif, A.,
Patel, C. J., & Members of the Exposomics Consortium. (2024). Decoding
the exposome: data science methodologies and implications in
exposome-wide association studies (ExWASs). *Exposome, 4*(1), osae001.
<https://doi.org/10.1093/exposome/osae001>

fenr. (n.d.). Bioconductor. Retrieved August 18, 2025, from
<https://www.bioconductor.org/packages/release/bioc/html/fenr.html>

Maitre, L., Guimbaud, J.-B., Warembourg, C., Güil-Oumrait, N., Petrone,
P. M., Chadeau-Hyam, M., Vrijheid, M., Basagaña, X., Gonzalez, J. R., &
Exposome Data Challenge Participant Consortium. (2022). State-of-the-art
methods for exposure-health studies: Results from the exposome data
challenge event. *Environment International, 168*(107422), 107422.
<https://doi.org/10.1016/j.envint.2022.107422>

Miller, G. W., & Banbury Exposomics Consortium. (2025). Integrating
exposomics into biomedicine. *Science, 388*(6745), 356–358.
<https://doi.org/10.1126/science.adr0544>

mixOmics. (n.d.). Bioconductor. Retrieved August 18, 2025, from
<https://www.bioconductor.org/packages/devel/bioc/html/mixOmics.html>

MOFA2. (n.d.). Bioconductor. Retrieved August 18, 2025, from
<https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html>

MultiAssay Special Interest Group. (2025, April 15).
*MultiAssayExperiment: The Integrative Bioconductor Container.*
<https://www.bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html>

nipalsMCIA. (n.d.). Bioconductor. Retrieved August 18, 2025, from
<https://www.bioconductor.org/packages/release/bioc/html/nipalsMCIA.html>

Regularized and Sparse Generalized Canonical Correlation Analysis for
Multiblock Data. (n.d.). Retrieved August 18, 2025, from
<https://rgcca-factory.github.io/RGCCA/>

Wild, C. P. (2005). Complementing the genome with an “exposome”: the
outstanding challenge of environmental exposure measurement in molecular
epidemiology. *Cancer Epidemiology, Biomarkers & Prevention, 14*(8),
1847–1850. <https://doi.org/10.1158/1055-9965.EPI-05-0456>
