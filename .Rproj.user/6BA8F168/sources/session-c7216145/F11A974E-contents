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
# inspiration:
# https://github.com/tidyomics
# https://stemangiola.github.io/tidybulk/articles/introduction.html

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
  "cd16_rna" = cd16_rna,
  "cd4_rna" = cd4_rna,
  "cd4_mirna"= cd4_mirna,
  "cd16_mirna" = cd16_mirna,
  "protein" = prot_abd_unnorm
)

expom <- expOmicSet(
  var_info = des,
  exposure = aw_dataset_cv2_filt,
  omics = omics_list,
  row_data = NULL)

colData(expom)
assays(expom)

exp_vars <- expom@metadata$var_info |> filter( category %in% c("allergen_results","chemical","indoor_pollution","Allergens","smoke_exposure")) |> pull(variable)

# --- Missingness Check --------
expom_1 <- expom |> 
  filter_missing(na_thresh = 20)

# --- Imputation --------
expom_2 <- expom_1 |> 
  impute_missing(
    exposure_impute_method = "median",
    omics_impute_method = "median")

# --- PCA Analysis --------
expom_3 <- expom_2 |> 
  pca_analysis()

# --- Data Normality Check ------------
expom_4 <- expom_3 |>
  check_normality()

# --- Transform Exposure Data ---------
expom_5 <- expom_4 |> 
  transform_exposure(transform_method = "best")

# --- Filter Non-Normal Exposure Data ---------
# expom_6 <- expom_5 |> 
#   filter_non_normal(p_thresh = 0.05)

# --- Subject Exposure Clustering --------------
expom_7 <- expom_5 |> 
  cluster_samples(exposure_cols = exp_vars,
                  clustering_approach = "diana")

# --- Exposure Correlation --------------------------
expom_8 <- expom_7 |>
  exposure_correlation()

# --- Exposure-Outcome Association -----------------
expom_9 <- expom_8 |> 
  perform_exwas(
    outcome = "pftfev1fvc_actual",
    exposures = colData(expom_8) |> 
      colnames() |> 
      (\(x) x[!x %in% c("pftfev1fvc_actual","age","black","female","income5")])(),
    confounders = c("age","black","female"))

# --- ExWAS Variable Selection ---------------------
# expom_10 <- expom_9 |>
#   exwas_select(
#     outcome = "pftfev1fvc_actual",
#     exposures = colData(expom_9) |>
#       colnames() |>
#       (\(x) x[!grepl("fev|fef|fvc|income5|black",x)])(),
#     confounders = c("age","black","female","income5"))

# --- Adjust Omics Data ---------------
# expom_11 <- expom_9 |>
#   adjust_assays(
#     outcome = "fev1fvc_category",
#     covariates = c("age","female", "black"),
#     minimum_counts = 10,
#     minimum_proportion = 0.7,
#     scaling_method = "none",
#     skip_identify_abundant = c("cd4_mirna","cd16_mirna","protein"))

# expom_11 <- expom_10 |>
#   adjust_assays(
#     outcome = "pftfev1fvc_actual",
#     covariates = c("age","female", "black"),
#     minimum_counts = 10,
#     minimum_proportion = 0.8,
#     scaling_method = "none",
#     skip_identify_abundant = c("cd4_mirna","cd16_mirna","protein"))

# --- Differential Analysis -------------
# a <- expom_9 |>
#   run_differential_abundance(
#     formula = ~ 0 + fev1fvc_category + age + sex + race,
#     minimum_counts = 1,
#     minimum_proportion = 0.1,
#     method = "deseq2",
#     contrasts = c("fev1fvc_categorySevere - fev1fvc_categoryUnobstructed",
#                   "fev1fvc_categoryModerate - fev1fvc_categoryUnobstructed"),
# 
#     scaling_method = "none")
# 
# a <- expom_9 |>
#   run_differential_abundance(
#     formula = ~ fev1fvc_category + age + sex + race,
#     minimum_counts = 1,
#     minimum_proportion = 0.1,
#     method = "limma_voom",
#     # contrasts = c("fev1fvc_categorySevere - fev1fvc_categoryUnobstructed",
#     #               "fev1fvc_categoryModerate - fev1fvc_categoryUnobstructed"),
#     contrasts = c("fev1fvc_categorySevere - (Intercept)",
#                   "fev1fvc_categoryModerate - (Intercept)"),
#     scaling_method = "none")

expom_10 <- expom_9 |>
  run_differential_abundance(
    formula = ~ pftfev1fvc_actual + age + sex + race,
    method = "limma_voom",
    contrasts = c("pftfev1fvc_actual - (Intercept)"),
    minimum_counts = 1,
    minimum_proportion = 0.7,
    scaling_method = "none")

expom_10@metadata[["differential_abundance"]] |>
  #filter(grepl("Severe",contrast)) |> 
  filter(adj.P.Val<0.05 & abs(logFC) > log2(1.5)) |> 
  #filter(adj.P.Val<0.05 ) |> 
  janitor::tabyl(assay_name) |> 
  arrange(desc(n))

# --- Feature:Exposure Correlation Analysis ----------
# expom_13 <- expom_12 |>
#   correlate_exposures_with_omics(
#     exposure_cols = exp_vars,
#     da_column = "adj.P.Val",
#     da_threshold = 0.05,
#     cor_pval_column = "FDR",
#     pval_cutoff = 0.1,
#     correlation_cutoff = 0.3
#   )

# expom_14 <- expom_13 |>
#   exposure_omic_association(
#     # exposures = colData(expom_12) |> 
#     #   colnames() |> 
#     #   (\(x) x[!x %in% c("pftfev1fvc_actual","age","race","sex")])(), 
#     exposures = exp_vars,
#     confounders = c("age","sex","race"), 
#     family = "gaussian", 
#     correction_method = "fdr"
#   )
# --- Sensitivity Analysis --------------
expom_11 <- run_sensitivity_analysis(
  expOmicSet = expom_10, 
  base_formula = ~ pftfev1fvc_actual + age + sex + race, 
  contrasts = c("pftfev1fvc_actual - (Intercept)"),
  methods = c("limma_voom"),
  scaling_methods = c("none"),
  min_counts_range = c(1, 5),
  min_proportion_range = c(0.1, 0.5),
  covariates_to_remove = c("age" , "sex" , "race")
)

# --- Multiomics Analysis -------------------
expom_12 <- expom_11 |> 
  multiomics_integration(method = "MCIA")

# --- Identify Relevant Factors -------------
expom_13 <- identify_relevant_factors(expOmicSet = expom_12, 
                                      outcome_var = "pftfev1fvc_actual", 
                                      categorical = FALSE,
                                      p_thresh = 0.1)

# --- Top Features By Factor Loadings ------
expom_14 <- expom_13 |> 
  extract_top_factor_features(factors = "V3", 
                              method = "percentile",
                              percentile = 0.95,
                              threshold = 0.3)

# --- Exposure-Omic Corrrelation ----------
expom_15 <- expom_14 |> 
  correlate_exposures_with_features(
    exposure_cols = c(exp_vars[!grepl("-",exp_vars)],"pftfev1fvc_actual"),
    batch_size = 1500)

# --- Functional Enrichment --------------
expom_16 <- expom_15 |> 
  multiomics_functional_enrichment(
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

