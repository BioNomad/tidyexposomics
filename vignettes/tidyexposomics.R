## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----pkg-overview, echo=FALSE,out.width="150%",fig.align='center',fig.cap="tidyexposomics pipeline overview. QC, association testing, integration, and enrichment steps on a MultiAssayExperiment."----
knitr::include_graphics("./overview.png")

## ----install, eval=FALSE------------------------------------------------------
# # install the package
# BiocManager::install("tidyexposomics")
# 
# # load the package
# library(tidyexposomics)

## ----command-str, echo=FALSE,fig.align='center',fig.cap= "Command naming conventions used throughout `tidyexposomics`. More complex pipelines begin with the `run_` prefix, visualizations with `plot_`, and data processing with `filter_`, `transform_`, `pivot_`, or `extract_` prefixes.",out.width="60%"----

knitr::include_graphics("./command_str.png")

## ----mae-overview,echo=FALSE,fig.align='center',fig.cap="Overview of the MultiAssayExperiment structure linking samples, assays, and metadata.",out.width="150%"----
knitr::include_graphics("./mae_rep.png")

## ----load-libraries,message=FALSE,warning=FALSE-------------------------------
# Load Libraries
library(tidyverse)
library(tidyexposomics)

## ----load-data----------------------------------------------------------------
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

## ----exp-vars-----------------------------------------------------------------
# Grab exposure variables
exp_vars <- tidyexposomics_example$annotated_cb |>
    filter(category %in% c(
        "aerosol",
        "main group molecular entity",
        "polyatomic entity"
    )) |>
    pull(variable) |>
    as.character()

## ----missing-bar,fig.height=3,fig.width=6,fig.cap="Count of features with missing data above a 0% missingness threshold by data layer. Exposure data have variables with missingness."----

# Plot the missingness summary
plot_missing(
    exposomicset = expom,
    plot_type = "summary",
    threshold = 0
)

## ----missing-bar-lollipop,fig.height=2,fig.width=4,fig.cap= "Percent missingness per exposure variable. Parity, `h_parity_None`, shows the highest missingness."----

# Plot missing variables withing exposure group
plot_missing(
    exposomicset = expom,
    plot_type = "lollipop",
    threshold = 0,
    layers = "Exposure"
)

## ----impute-missing-----------------------------------------------------------

# Impute missing values
expom <- run_impute_missing(
    exposomicset = expom,
    exposure_impute_method = "missforest",
    exposure_cols = exp_vars
)

## ----filter-omics-------------------------------------------------------------

# filter omics layers by variance and expression
# Methylation filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "variance",
    assays = "Methylation",
    assay_name = 1,
    min_var = 0.05
)

# Gene expression filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "expression",
    assays = "Gene Expression",
    assay_name = 1,
    min_value = 1,
    min_prop = 0.3
)

## ----normality-check----------------------------------------------------------

# Check variable normality
expom <- run_normality_check(
    exposomicset = expom,
    action = "add"
)

## ----trasform-vars------------------------------------------------------------

# Transform variables
expom <- transform_exposure(
    exposomicset = expom,
    transform_method = "boxcox_best",
    exposure_cols = exp_vars
)

## ----norm-plot,fig.width= 4,fig.height= 4, fig.cap= "Normality status of numeric exposure variables after Box-Cox transformation."----

# Examine normality summary
plot_normality_summary(
    exposomicset = expom,
    transformed = TRUE
)

## ----pca----------------------------------------------------------------------

# Perform principal component analysis
expom <- run_pca(
    exposomicset = expom,
    log_trans_exp = TRUE,
    log_trans_omics = TRUE,
    action = "add"
)

## ----pca-plot, fig.align='center', fig.width= 9, fig.height= 7,fig.cap= "PCA of sample and feature space with sample outlier detection."----

# Plot the PCA plot of sample and feature space
plot_pca(exposomicset = expom)

## ----filt-outliers------------------------------------------------------------

# Filter out sample outliers
expom <- filter_sample_outliers(
    exposomicset = expom,
    outliers = c("s1231")
)

## ----corr-pcs-----------------------------------------------------------------

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

## ----plot-pc-corr,fig.width= 7,fig.height= 2.2,fig.cap= "Correlation heatmap of exposures versus principal components. Child lead levels (`hs_pb_c_Log2`) and maternal BPA levels (`hs_bpa_madj_Log2`) are associated with the most principal components."----


# Plot the correlation tile plot
plot_correlation_tile(
    exposomicset = expom,
    feature_type = "pcs",
    pval_cutoff = 0.05
)

## ----exp-sum, fig.width=6, fig.height=4---------------------------------------

# Summarize exposure data
run_summarize_exposures(
    exposomicset = expom,
    action = "get",
    exposure_cols = exp_vars
) |>
    head()

## ----plot-aerosol, fig.width=4, fig.height=3.5, fig.cap="Distribution of aerosol exposures by sex."----

# Plot aerosol exposure distributions by sex
plot_exposures(
    exposomicset = expom,
    group_by = "e3_sex_None",
    exposure_cat = "aerosol",
    plot_type = "boxplot",
    ylab = "Values",
    title = "Aerosol Exposure by Sex"
)

## ----sample-clusters----------------------------------------------------------

# Sample clustering
expom <- run_cluster_samples(
    exposomicset = expom,
    exposure_cols = exp_vars,
    clustering_approach = "dynamic",
    action = "add"
)

## ----plot-sample-clusters, fig.height=4, fig.width=12, fig.cap="Sample clustering heatmap using exposure profiles (z-scored). Clusters appear mostly driven by aerosol exposure during pregnancy."----

# Plot the sample clusters
plot_sample_clusters(
    exposomicset = expom,
    exposure_cols = exp_vars
)

## ----correlate-exposures------------------------------------------------------

# Run correlation analysis
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "exposures",
    action = "add",
    exposure_cols = exp_vars,
    correlation_cutoff = 0.3
)

## ----exposure-circos-corr, fig.height=10,fig.width=10,fig.align='center',fig.cap="Circos view of exposure-exposure correlations (threshold 0.3)."----

# Plot exposure correlation circos plot
plot_circos_correlation(
    exposomicset = expom,
    feature_type = "exposures",
    corr_threshold = 0.3,
    exposure_cols = exp_vars
)

## ----assoc-exposures----------------------------------------------------------

# Perform ExWAS
expom <- run_association(
    exposomicset = expom,
    source = "exposures",
    outcome = "hs_asthma",
    feature_set = exp_vars,
    action = "add",
    family = "binomial"
)

## ----plot-exposure-assoc,fig.height=3, fig.width=5, fig.align='center',fig.cap="ExWAS associations of exposures with asthma status. No exposures are significantly associated (P < 0.05) with asthma status."----

# Plot the association forest plot
plot_association(
    exposomicset = expom,
    source = "exposures",
    terms = exp_vars,
    filter_thresh = 0.05,
    filter_col = "p.value",
    r2_col = "r2"
)

## ----associate-top-omics,fig.height=3,fig.width=4.5,fig.align='center'--------

# Perform ExWAS
expom <- run_association(
    exposomicset = expom,
    outcome = "hs_asthma",
    source = "omics",
    top_n = 500,
    action = "add",
    family = "binomial"
)

## ----plot-manhattan,fig.height=7, fig.width=5, fig.align='center',fig.cap="Manhattan plot of omics-wide associations with asthma status."----

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

## ----diff-abundance,results='hide'--------------------------------------------

# Run differential abundance analysis
expom <- run_differential_abundance(
    exposomicset = expom,
    formula = ~hs_asthma,
    method = "limma_trend",
    scaling_method = "none",
    action = "add"
)

## ----volcano-plot,fig.height=3.5, fig.width=6, fig.align='center',fig.cap="Volcano plot of differentially abundant features across omics layers."----

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

## ----exp-omics-assoc, fig.align='center'--------------------------------------
# Run association testing between every exposure and omics feature
expom <- run_exposure_omics_association(
    exposomicset = expom,
    exposures = exp_vars,
    action = "add"
)

## ----plot-exp-omics-assoc, fig.align='center',fig.width=7,fig.height=2.5,fig.cap="Barplot of the number of omics features associated with exposures."----

# Plot the number of exposure-omics associations
plot_exposure_omics_association(
  exposomicset = expom,
  plot_type = "category",
  pval_col = "p.value",
  pval_thresh = 0.05)

## ----var-map------------------------------------------------------------------
# Extract omics features associated with aerosols
var_map <- extract_results(
  exposomicset = expom,
  result = "association"
) |> 
  pluck("exposure_omics",
        "results_df") |> 
  filter(p.value<0.05) |>
  filter(category == "aerosol") |> 
  dplyr::select(
    exp_name = exp_name,
    variable = feature
    
  )

## ----enrichment---------------------------------------------------------------
# Run enrichment analysis on factor features correlated with exposures
expom <- run_enrichment(
    exposomicset = expom,
    variable_map = var_map,
    feature_type = "omics",
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

## ----enrichment-summary,fig.width=14,fig.height=5, fig.cap="Summary of enriched GO terms grouped by overlap and exposure category."----

# Plot the summary diagram
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics",
    plot_type = "summary"
)

## ----enrichment-dotplot,fig.height=11, fig.width=11, fig.cap="Dotplot of top enriched GO terms by omics layer and exposure category."----

# Plot a dotplot of terms
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics",
    plot_type = "dotplot",
    top_n = 15,
    add_top_genes = TRUE,
    top_n_genes = 5
)

## ----enrichment-term-network, fig.width=6, fig.height=5, fig.cap="Network of enriched GO terms connected by shared genes."----

# Plot the term network plot
# Setting a seed so that the plot layout is consistent
set.seed(42)
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics",
    plot_type = "network",
    label_top_n = 3
)

## ----enrichment-heatmap, warning=FALSE,fig.height=3, fig.width=8, fig.cap="Heatmap of genes driving enriched GO terms (Group 2) with log2 fold-change overlay."----

# Plot a heatmap of genes and corresponding GO terms
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics",
    go_groups = "Group 2",
    plot_type = "heatmap",
    heatmap_fill = TRUE,
    feature_col = "feature_clean"
)

## ----enrichment-cnetplot, fig.width=6, fig.height=6, fig.cap="Cnet plot linking enriched terms to contributing genes (Group 1)."----

# Plot the gene-term network
# Setting a seed so that the plot layout is consistent
set.seed(42)
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics",
    go_groups = "Group 2",
    plot_type = "cnet",
    feature_col = "feature_clean"
)

## ----pipeline-summary---------------------------------------------------------

# Run the pipeline summary
expom |>
    run_pipeline_summary(console_print = TRUE, include_notes = TRUE)

## ----save-results-------------------------------------------------------------

# Save results
extract_results_excel(
    exposomicset = expom,
    file = tempfile(),
    result_types = "all"
)

## ----session_info,collapse = TRUE---------------------------------------------

sessionInfo()

