# --- Testing tidyexposomics --------------
library(tidyverse)
library(tidyomics)
library(readxl)
library(clusterProfiler)
invisible(lapply(
  list.files(path = "~/jhu/projects/tidyexposomics/R/",
             pattern="*.R",full.names = TRUE),
  source))

# --- Testing tidyexposomics --------------

# --- Load Data ---------------------------
load("./data/expom.RData")

# TODO: check s2809 - manually coding as male based on gene exp of Y chr genes
aw_dataset_cv2_filt <- aw_dataset_cv2_filt |> 
  mutate(female=case_when(
    aw_id == "s2809" ~ "male",
    .default=female
  )) |> 
  mutate(race=case_when(
    black=="black" ~ "black",
    black=="non-black" ~ "non_black",
  )) |> 
  mutate(fev1fvc_category=case_when(
    fev1fvc_category == "Normal" ~ "Unobstructed",
    fev1fvc_category == "Moderate" ~ "Moderate",
    fev1fvc_category == "Severe" ~ "Severe"
  )) |>
  mutate(fev1fvc_category=factor(fev1fvc_category,levels=c("Unobstructed","Moderate","Severe"))) |> 
  mutate(age_scaled=as.numeric(scale(age))) |> 
  mutate(sex=female) |> 
  mutate(sex=factor(sex,levels=c("female","male"))) |>
  mutate(race=factor(race,levels=c("non_black","black")))

# protein is low variance, rescaling
prot_abd_log <- prot_abd |> 
  (\(x) log2(x+1))() |> 
  as.data.frame()

prot_abd_unnorm <- prot_abd |> 
  (\(x) 2^x)()

omics_list <- list(
  "CD16+ Monocyte RNA" = cd16_rna,
  "CD4+ T-cell RNA" = cd4_rna,
  "CD4+ T-cell miRNA"= cd4_mirna,
  "CD16+ Monocyte miRNA" = cd16_mirna,
  "Serum Proteomics" = prot_abd_unnorm
)

expom <- expOmicSet(
  var_info = des,
  exposure = aw_dataset_cv2_filt,
  omics = omics_list,
  row_data = NULL)

# option if you have one omic
# expom <- expOmicSet(
#   var_info = des,
#   exposure = aw_dataset_cv2_filt,
#   omics = list("cd4_rna"=omics_list$cd4_rna),
#   row_data = NULL)

colData(expom)
assays(expom)

exp_vars <- expom@metadata$var_info |> 
  filter( category %in% 
            c("allergen_results",
              "chemical",
              "indoor_pollution",
              "Allergens",
              "smoke_exposure")) |> 
  pull(variable)

exp_vars <- c(exp_vars[exp_vars %in% colnames(colData(expom))])

# --- Missingness Check --------
expom <- expom |> 
  filter_missing(na_thresh = 20)

expom |> 
  plot_missing_summary(threshold = 20)

# --- Imputation --------
expom <- expom |> 
  impute_missing(
    exposure_impute_method = "median",
    omics_impute_method = "median")

# --- PCA Analysis --------
expom <- expom |> 
  pca_analysis()

expom |> 
  plot_pca()

# a <- expom |> 
#   remove_sample_outliers()
# 
# a <- a |> 
#   pca_analysis()
# 
# a |> 
#   plot_pca()

# --- Data Normality Check ------------
expom <- expom |>
  check_normality()

# --- Transform Exposure Data ---------
expom <- expom |> 
  transform_exposure(transform_method = "best")

expom |> 
  plot_normality_sum()

# --- Filter Non-Normal Exposure Data ---------
# expom_6 <- expom_5 |> 
#   filter_non_normal(p_thresh = 0.05)

# --- Subject Exposure Clustering --------------
expom <- expom |> 
  cluster_samples(exposure_cols = exp_vars,
                  clustering_approach = "diana")

expom |> 
  plot_sample_clusters(
    cols_of_interest = c(exp_vars,"age","pftfev1fvc_actual")) 

# --- Exposure Correlation --------------------------
expom <- expom |>
  exposure_correlation(
    exposure_cols = c(exp_vars,"age","pftfev1fvc_actual"))

expom |> 
  plot_exposure_corr()

# --- Exposure-Outcome Association -----------------
a <- expom@metadata$var_info |> filter(category != "spirometry") |> pull(variable)
expom <- expom |> 
  perform_exwas(
    outcome = "pftfev1fvc_actual",
    exposures = a[!a %in% c("pftfev1fvc_actual",
                            "age",
                            "sex",
                            "race")],
    # exposures = colData(expom) |> 
    #   colnames() |> 
    #   (\(x) x[!x %in% 
    #             c("pftfev1fvc_actual",
    #               "age",
    #               "black",
    #               "female",
    #               "income5")])(),
    covariates = c("age","sex","race"))

expom |> 
  plot_exp_outcome_association()

# --- Differential Analysis -------------
expom <- expom |>
  run_differential_abundance(
    formula = ~ pftfev1fvc_actual + age + sex + race,
    method = "limma_voom",
    contrasts = c("pftfev1fvc_actual - (Intercept)"),
    minimum_counts = 1,
    minimum_proportion = 0.7,
    scaling_method = "none")

expom |> 
  plot_volcano()

expom@metadata[["differential_abundance"]] |>
  filter(adj.P.Val<0.05 & abs(logFC) > log2(1.5)) |> 
  janitor::tabyl(assay_name) |> 
  arrange(desc(n))

# --- Sensitivity Analysis --------------
expom <- expom |> 
  run_sensitivity_analysis(
  base_formula = ~ pftfev1fvc_actual + age + sex + race, 
  contrasts = c("pftfev1fvc_actual - (Intercept)"),
  methods = c("limma_voom"),
  scaling_methods = c("none"),
  min_counts_range = c(1,5),
  min_proportion_range = c(0.1,0.5),
  covariates_to_remove = c("age" , "sex" , "race")
)

expom |> 
  plot_sensitivity_summary()

# --- Multiomics Analysis -------------------
expom <- expom |> 
  multiomics_integration(method = "MCIA")

expom |> 
  plot_factor_summary()

expom |> 
  plot_top_factor_features()

a <- expom |> 
  multiomics_integration(method = "MOFA")

# --- Identify Relevant Factors -------------
expom <- expom |> 
  identify_relevant_factors(outcome_var = "pftfev1fvc_actual", 
                                      categorical = FALSE,
                                      p_thresh = 0.1)

a <- a |> 
  identify_relevant_factors(outcome_var = "pftfev1fvc_actual", 
                            categorical = FALSE,
                            p_thresh = 0.1)

# --- Top Features By Factor Loadings ------
expom <- expom |> 
  extract_top_factor_features(factors = "V3", 
                              method = "percentile",
                              percentile = 0.95,
                              threshold = 0.3)

a <- a |> 
  extract_top_factor_features(factors = "Factor1", 
                              method = "percentile",
                              percentile = 0.95,
                              threshold = 0.3)

# --- Feature:Exposure Correlation Analysis ----------
expom <- expom |> 
  correlate_exposures_with_degs(exposure_cols = exp_vars)

expom <- expom |> 
  correlate_exposures_with_factors(exposure_cols = exp_vars)

# --- Functional Enrichment --------------
expom <- expom |> 
  .exposure_category_functional_enrichment(
    cor_df = "degs",
    mirna_assays = c("cd4_mirna", "cd16_mirna"), 
    uniprot_assays = c("protein"))


# --- TODO ---------------------------
#| - Ensure that each type of enrichment algorithm enriches each omic separately
#| - Determine what are the shared terms and what are different per omic
#| - Determine what are the shared terms and what are different per exposure category
#| - Add in session info to the expOmicSet object
#| - Create plotting functions for each step
#| - Create functions that can export the results to an excel file or csv
#| - Change the correlation function to apply FDR adjustment at the end and not per chunk
#| - For each of the enrichment terms, cluster them, what groups are enriched in each omic and for each exposure category, what is different and what is the same?
#| Code to use:upset(fromList( map(b@metadata$exposure_category_enrichment,~ .x |> filter(p.adjust<0.05) |> pull(Description))),nintersects = 100,nsets = 10) 
#| 
# --- Pathway Analysis ----------------
expom_19 <- expom_18 |> 
  multiomics_pathway_analysis(
    omic = "cd4_rna",
    factors = 3)

# --- Exposure-Omic Corrrelation ----------
expom_19 <- expom_18 |> 
  correlate_exposures_with_features(exposure_cols = c(exp_vars[!grepl("-",exp_vars)],"pftfev1fvc_actual"))

# --- Factor-Exposure Correlation -------
# expom_16 <- expom_15 |> 
#   correlate_factors_with_exposures(exposures = c(exp_vars[!grepl("-",exp_vars)],"pftfev1fvc_actual"))

# --- Visualization ----------------------

# --- Overflow -----------
# --- Volcano Plot Per Omic ------------
expom_13@metadata$differential_abundance |> 
  mutate(assay_name=case_when(
    assay_name == "cd16_rna" ~ "CD16+ Monocyte RNA",
    assay_name == "cd4_rna" ~ "CD4+ T-cell RNA",
    assay_name == "cd4_mirna" ~ "CD4+ T-cell miRNA",
    assay_name == "cd16_mirna" ~ "CD16+ Monocyte miRNA",
    assay_name == "protein" ~ "Serum Protein"
  )) |> 
  mutate(direction=case_when(
    ((logFC > 1.5)  & (adj.P.Val < 0.05)) ~ "Up-Regulated",
    ((logFC < -1.5) & (adj.P.Val < 0.05)) ~ "Down-Regulated",
    .default = "No Change")) |>
  ggplot(aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = direction))+
  geom_point(alpha=.7)+
  theme_pubr(legend="bottom")+
  scale_color_manual(values=c(
    "No Change" = "grey65",
    "Up-Regulated" = "magenta4",
    "Down-Regulated" = "midnightblue"))+
  labs(
    x = expression("Log"[2]*" Fold Change"),
    y = expression("-Log"[10]*" P"),
    color = "")+
  facet_wrap(~assay_name,scales = "free",nrow = 1)+
  theme(strip.text = element_text(face="bold.italic"),
        axis.text = element_text(face="italic"))

# --- Enrichment Analysis ------------------------

# read in protein mapped to gene
pmap <- read_excel("../exposome/data/proteomics.xlsx") |> 
  (\(x) x[,1:2])()

allergen_da <- expom_13@metadata$omics_exposure_correlation |>  
  filter(grepl("allergen",category)) |> 
  left_join(pmap, 
            by=c("omics"="Protein.Group")) |> 
  mutate(omics=case_when(
    is.na(Genes) ~ omics,
    !is.na(Genes) ~ Genes)) 

chemical_da <- expom_13@metadata$omics_exposure_correlation |>  
  filter(grepl("chemical",category)) |> 
  left_join(pmap, 
            by=c("omics"="Protein.Group")) |> 
  mutate(omics=case_when(
    is.na(Genes) ~ omics,
    !is.na(Genes) ~ Genes))

air_pollution_da <- expom_13@metadata$omics_exposure_correlation |>  
  filter(grepl("indoor_pollution|smoke_exposure",category)) |> 
  left_join(pmap, 
            by=c("omics"="Protein.Group")) |> 
  mutate(omics=case_when(
    is.na(Genes) ~ omics,
    !is.na(Genes) ~ Genes))

exp_genes <- list(
  "allergen" = allergen_da |> 
    filter(FDR<0.05) |> 
    pull(omics) |> 
    unique(),
  "chemical" = chemical_da |> 
    filter(FDR<0.05) |>  
    pull(omics) |> 
    unique(),
  "air_pollution" = air_pollution_da |> 
    filter(FDR<0.05) |> 
    pull(omics) |> 
    unique()
)

# apply clusterProfilter to each list element
enrich_res <- lapply(exp_genes, function(genes) {
  res <- enrichGO(
    gene = genes,
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
})

enrich_res_df <- map2(enrich_res, names(enrich_res), function(res, name) {
  res@result |> 
    mutate(exposure = name) 
}) |> 
  bind_rows()

enrich_res_df |> 
  filter(grepl("oxidative|stress|antigen",Description)) |> 
  mutate(exposure = case_when(
    exposure == "allergen" ~ "Allergen",
    exposure == "chemical" ~ "Chemical",
    exposure == "air_pollution" ~ "Air Pollution"
  )) |> 
  group_by(exposure) |> 
  arrange(p.adjust) |> 
  slice_head(n=10) |> 
  ungroup() |> 
  ggplot(aes(
    x = exposure,
    y = reorder(Description, -log10(p.adjust)),
    size = Count,
    color = -log10(p.adjust)))+
  geom_point()+
  theme_pubr(legend="right")+
  rotate_x_text(angle=45)+
  scale_color_gradient(low="thistle",high = "midnightblue")+
  labs(
    x = "",
    y = "",
    title = paste("Enriched Gene Ontology Terms", "By Exposure Type",sep="\n"),
    color = expression("-Log"[10]*" P"))+
  theme(plot.title = element_text(face="bold.italic"))

expom_13@metadata$omics_exposure_correlation |>  
  left_join(pmap, 
            by=c("omics"="Protein.Group")) |> 
  mutate(omics=case_when(
    is.na(Genes) ~ omics,
    !is.na(Genes) ~ Genes)) |> 
  mutate(category=gsub("smoke_exposure","indoor_pollution",category)) |> 
  janitor::tabyl(omics,category) |> 
  pivot_longer(
    cols = c("allergen_results","chemical","indoor_pollution"),
    names_to = "exposure",
    values_to = "n") |> 
  group_by(omics) |> 
  mutate(total=sum(n)) |> 
  ungroup() |> 
  arrange(desc(total)) |> 
  filter(omics %in% (expom_13@metadata$omics_exposure_correlation |>  
                       left_join(pmap, 
                                 by=c("omics"="Protein.Group")) |> 
                       mutate(omics=case_when(
                         is.na(Genes) ~ omics,
                         !is.na(Genes) ~ Genes)) |> janitor::tabyl(omics) |> arrange(desc(n)) |> slice_head(n=50) |> pull(omics))) |> 
  mutate(exposure = case_when(
    exposure == "allergen_results" ~ "Allergen",
    exposure == "chemical" ~ "Chemical",
    exposure == "indoor_pollution" ~ "Air Pollution"
  )) |> 
  ggplot(aes(
    x = n,
    y = fct_reorder(omics, total),
    fill=exposure
  )) +
  geom_bar(stat = "identity",alpha = 0.8) +
  theme_pubr(legend = "right") +
  scale_color_cosmic(guide="none")+
  scale_fill_cosmic()+
  labs(x = paste("No. of Exposure", "Associations",sep="\n"),
       y = "",fill="",
       title = paste("No. of Exposure Associations By", "Molecular Feature",sep="\n"))+
  theme(plot.title = element_text(face="bold.italic"))

# --- Correlation Network Analysis ------------
# not sure about this one
# expom_14 <- expom_13 |> 
#   run_correlation_network(
#     assay_to_use = 1,
#     merge_threshold = 0.75)

