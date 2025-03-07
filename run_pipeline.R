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

expom <- create_expomicset(
  var_info = des,
  exposure = aw_dataset_cv2_filt,
  omics = omics_list,
  row_data = fdata)

expom_qc <- expom |> 
  # Filter out variables with too many missing values
  filter_missing(na_thresh = 20) |>
  
  # Impute missing values
  impute_missing(
    exposure_impute_method = "median",
    omics_impute_method = "median") |> 
  
  # Principal component analysis
  pca_analysis() |> 
  
  # Check Variable Normality
  check_normality() 
  
  # Transform Variables 
  transform_exposure(transform_method = "best") 

expom_sample_exp <- expom_qc |> 
  # Sample Clustering
  cluster_samples(exposure_cols = exp_vars,
                  clustering_approach = "diana") |> 
  
  # Perform ExWAS Analysis
  perform_exwas(
    outcome = "fev_height",
    exposures = all_vars[
      !all_vars %in% c("fev_height",
                       "age",
                       "sex",
                       "race",
                       "smoking")],
    covariates = c("age","sex","race","smoking")) 

expom_da <- expom_sample_exp |> 
  # Perform Differential Abundance Analysis
  run_differential_abundance(
    formula = ~ pftfev1fvc_actual + age + sex + race,
    method = "limma_voom",
    contrasts = c("pftfev1fvc_actual - (Intercept)"),
    minimum_counts = 1,
    minimum_proportion = 0.1,
    scaling_method = "none") 

expom_da <- expom_da |> 
  # Perform Sensitivity Analysis
  run_sensitivity_analysis(
    base_formula = ~ pftfev1fvc_actual + age + sex + race, 
    contrasts = c("pftfev1fvc_actual - (Intercept)"),
    methods = c("limma_voom"),
    scaling_methods = c("none"),
    min_counts_range = c(1,5),
    min_proportion_range = c(0.1,0.3),
    covariates_to_remove = c("age" , "sex" , "race")) 

expom_multi <- expom_da |> 
  # Perform Multi-Omics Integration
  multiomics_integration(method = "MCIA") |> 
  
  # Identify factors that correlate with the outcome
  identify_relevant_factors(outcome_var = "pftfev1fvc_actual", 
                            categorical = FALSE,
                            p_thresh = 0.1) |> 
  extract_top_factor_features(factors = "V3", 
                              method = "percentile",
                              percentile = 0.95,
                              threshold = 0.3) 

expom_multi <- expom_multi |> 
  
  # Identify DEGs that correlate with exposures
  correlate_exposures_with_degs(exposure_cols = exp_vars) |> 
  
  # IDentify DEGs that correlate with factors
  correlate_exposures_with_factors(exposure_cols = exp_vars) 


expom_enrich <- expom_multi |> 
  # Perform Functional Enrichment Analysis
  run_functional_enrichment(
    geneset = "deg_exp_cor",
    feature_col = "gene",
    mirna_assays = c("CD16+ Monocyte miRNA","CD4+ T-cell miRNA"),
    pval_threshold = 0.05,
    logfc_threshold = log2(1.5))

# --- Test Functionality -------------------
expom_enrich |> pivot_sample()