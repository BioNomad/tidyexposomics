# --- Load Data ---------------------------
load("./data/expom.RData")
load("./data/exp_vars.RData")
load("./data/all_vars.RData")
invisible(lapply(
  list.files(path = "~/jhu/projects/tidyexposomics/R/",
             pattern="*.R",full.names = TRUE),
  source))

all_vars <- des |> filter(category!="spirometry") |> pull(variable)

exp_vars <- des |> filter(category %in% c(
  "urin_metal",
  "serum_metal",
  "indoor_pollution",
  "allergen_results",
  "chemical"
)) |> 
  pull(variable)

omics_list <- list(
  "CD16+ Monocyte RNA" = cd16_rna,
  "CD4+ T-cell RNA" = cd4_rna,
  "CD4+ T-cell miRNA"= cd4_mirna,
  "CD16+ Monocyte miRNA" = cd16_mirna,
  "Serum Proteomics" = prot_abd_unnorm,
  "Serum Adductomics" = adduct
)

fdata <- list(
  "CD16+ Monocyte RNA" = cd16_rna_fdata,
  "CD4+ T-cell RNA" = cd4_rna_fdata,
  "CD4+ T-cell miRNA"= cd4_mirna_fdata,
  "CD16+ Monocyte miRNA" = cd16_mirna_fdata,
  "Serum Proteomics" = prot_fdata,
  "Serum Adductomics" = adduct_fdata
)

# --- Create ExpomicSet -------------------
expom <- create_expomicset(
  var_info = des,
  exposure = aw_dataset_cv2_filt,
  omics = omics_list,
  row_data = fdata)

# --- Quality Control --------------------
expom_qc <- expom |> 
  # Filter out variables with too many missing values
  filter_missing(na_thresh = 20) |>
  
  # Impute missing values
  run_impute_missing(
    exposure_impute_method = "median",
    omics_impute_method = "median") |> 
  
  # Principal component analysis
  run_pca_analysis() |> 
  
  # Check Variable Normality
  run_normality_check() 
  
  # Transform Variables 
expom_qc <- expom_qc |> 
  transform_exposure(transform_method = "best",
                     cols_of_interest = exp_vars) 

# --- Sample Clustering and ExWAS Analysis ----
expom_sample_exp <- expom_qc |> 
  # Sample Clustering
  run_cluster_samples(exposure_cols = exp_vars,
                  clustering_approach = "dynamic") |> 
  
  # Perform ExWAS Analysis
  associate_exposure_outcome(
    outcome = "fev_height",
    exposures = all_vars[
      !all_vars %in% c("fev_height",
                       "age",
                       "sex",
                       "race")],
    covariates = c("age","sex","race")) 

# Plot Sample Clustering Results
expom_sample_exp |> plot_sample_clusters(cols_of_interest = exp_vars)

# Plot ExWAS Results
expom_sample_exp |> plot_associate_exposure_outcome(filter_col = "p.value",filter_thresh = 0.1)

# --- Differential Abundance Analysis ----
expom_da <- expom_sample_exp |> 
  # Perform Differential Abundance Analysis
  run_differential_abundance(
    formula = ~ fev_height + age + sex + race,
    method = "limma_voom",
    contrasts = c("fev_height - (Intercept)"),
    minimum_counts = 1,
    minimum_proportion = 0.1,
    scaling_method = "none") 

# Plot Differential Abundance Results
expom_da |> plot_volcano(logFC_thresh = log2(2),pval_thresh = 0.05)

# --- Sensitivity Analysis --------------------
expom_da <- expom_da |> 
  # Perform Sensitivity Analysis
  run_sensitivity_analysis(
    base_formula = ~ fev_height + age + sex + race, 
    contrasts = c("fev_height - (Intercept)"),
    methods = c("limma_voom"),
    scaling_methods = c("none"),
    min_counts_range = c(1,5),
    min_proportion_range = c(0.1,0.3),
    covariates_to_remove = c("age" , "sex" , "race")) 

# Plot sensitivity analysis results
expom_da |> plot_sensitivity_summary(stability_score_thresh = 12)

# --- Multi-Omics Integration --------------------
expom_multi <- expom_da |> 
  # Perform Multi-Omics Integration
  run_multiomics_integration(method = "MCIA") |> 
  
  # Identify factors that correlate with the outcome
  associate_factor_outcome(outcome_var = "fev_height", 
                            categorical = FALSE,
                            p_thresh = 1) 

# Extract top features that contribute to a factor
expom_multi <- expom_multi |> 
  extract_top_factor_features(factors = "V6", 
                              method = "percentile",
                              percentile = 0.95,
                              threshold = 0.3) 

# Plot multi-omics factor summary
expom_multi |> plot_factor_summary()

# --- Exposure-Omic Correlation Analysis --------------------
expom_multi <- expom_multi |> 
  
  # Identify DEGs that correlate with exposures
  correlate_exposures_degs(exposure_cols = exp_vars,
                                robust = TRUE,
                                score_thresh = 12) |> 
  
  # Identify DEGs that correlate with factors
  correlate_exposures_factors(exposure_cols = exp_vars) 

# Plot Exposure-Omic Correlation Results
expom_multi |> plot_bar_correlate_summary()

# Plot Shared Feature Correlations Between Exposures
expom_multi |> plot_circos_exposure_shared_features(geneset = "degs")

# --- Functional Enrichment Analysis --------------------
expom_enrich <- expom_multi |> 
  # Perform Functional Enrichment Analysis
  run_enrichment(
    geneset = "deg_exp_cor",
    feature_col = "gene",
    mirna_assays = c("CD16+ Monocyte miRNA","CD4+ T-cell miRNA"),
    pval_threshold = 0.05,
    logfc_threshold = log2(1.5))

expom_enrich <- expom_enrich |> 
  # Identify GO Groups that correlate with the outcome
  associate_go_outcome(
    geneset = "deg_exp_cor",
    outcome = "fev_height",
    mirna_assays = c("CD4+ T-cell miRNA","CD16+ Monocyte miRNA"),
    covariates = c("age","sex","race")
  )

# Plot GO Group Eigengene, Outcome Association
expom_enrich |> plot_associate_go_outcome(direction_filter = "down",filter_thresh = 0.15)

# Plot Functional Enrichment Results
expom_enrich |> plot_dotplot_enrichment(geneset = "deg_exp_cor",go_groups=c("3","8","28","6","5","20","11","29","18","16","17"))

# --- Test Functionality -------------------
expom_enrich |> pivot_sample()
expom_enrich |> pivot_feature()