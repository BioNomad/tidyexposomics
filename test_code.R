## tidyexposomics_vignette.R

# Load necessary libraries
library(tidyverse)
#library(tidyexposomics)

# Load data and source functions
load("./data/expom.RData")
invisible(lapply(
  list.files(path = "~/jhu/projects/tidyexposomics/R/", pattern = "*.R", full.names = TRUE),
  source))

# Create omics list and feature data
omics_list <- list(
  "Gene Expression" = exp,
  "Metabolomics" = met,
  "Proteomics" = prot
)

fdata <- list(
  "Gene Expression" = exp_fdata,
  "Metabolomics" = met_fdata,
  "Proteomics" = prot_fdata
)

# Create expomicset object
expom <- create_expomicset(
  var_info = des,
  exposure = meta,
  omics = omics_list,
  row_data = fdata)

exp_vars <- des |>
  filter(grepl("exposure|chemical", category, ignore.case = TRUE)) |>
  pull(variable) |>
  as.character()

# Run full tidyexposomics pipeline ----------------------------------------

# Quality control
expom <- expom |>
  filter_missing(na_thresh = 1) |>
  run_impute_missing(exposure_impute_method = "median", omics_impute_method = "median") |>
  run_pca(action = "add") |>
  filter_sample_outliers(outliers = c("s411","s378", "s588", "s764", "s857", "s918", "s936")) |>
  run_normality_check(action = "add") |>
  transform_exposure(transform_method = "boxcox_best", exposure_cols = exp_vars)

# Plot diagnostics
expom |> plot_missing_summary(threshold = 1)
expom |> plot_pca()
expom |> plot_normality_summary(transformed = TRUE)

# Summary and visualization
expom |> run_summarize_exposures(action = "get") |> head()
expom |> plot_exposures(panel_sizes = c(3,4), exposure_cols = exp_vars, plot_type = "boxplot", ylab = "Values")

# Clustering and correlation
expom <- expom |> run_cluster_samples(exposure_cols = exp_vars, clustering_approach = "dynamic", action = "add")
expom |> plot_sample_clusters(exposure_cols = exp_vars)
expom <- expom |> run_correlation(feature_type = "exposures", action = "add", correlation_cutoff = 0.3)
expom |> plot_circos_correlation(feature_type = "exposures", corr_threshold = 0.3, exposure_cols = exp_vars)

# Exposure association
expom <- expom |> run_association(source = "exposures", outcome = "hs_asthma",
                                  feature_set = exp_vars[!exp_vars %in% c("hs_asthma","hs_child_age","e3_sex")],
                                  covariates = c("hs_child_age", "e3_sex"), action = "add", family = "binomial")
expom |> plot_association(subtitle = "Covariates: Age, Sex", source = "exposures", filter_thresh = 0.05)

# Omics association
expom <- expom |> run_association(outcome = "hs_asthma", source = "omics",
                                  covariates = c("hs_child_age", "e3_sex"), top_n = 500, action = "add", family = "binomial")
expom |> plot_manhattan(min_per_cat = 0, vars_to_label = c("TC17001579.hg.1","TC02003949.hg.1","HMGN2P28"), panel_sizes = c(1,4,2,1,3))

# Exposome scores
scores <- c("median", "pca", "irt", "quantile", "var")
for (score in scores) {
  score_col <- paste0("exposome_", score, "_score")
  expom <- expom |> run_exposome_score(exposure_cols = exp_vars, score_type = score, score_column_name = score_col)
}
expom <- expom |> run_association(outcome = "hs_asthma", source = "exposures",
                                  feature_set = paste0("exposome_", scores, "_score"),
                                  covariates = c("hs_child_age", "e3_sex"), action = "add", family = "binomial")
expom |> plot_association(subtitle = "Covariates: Age, Sex", source = "exposures",
                          terms = paste0("exposome_", scores, "_score"), filter_thresh = 1)

# Differential abundance and sensitivity
expom <- expom |> run_differential_abundance(
  formula = ~ hs_asthma + hs_child_age + e3_sex,
  method = "limma_voom", minimum_counts = 1,
  minimum_proportion = 0.1, scaling_method = "none", action = "add")
expom |> plot_volcano(top_n_label = 3, logFC_thresh = log2(1), pval_thresh = 0.05, nrow = 1)
expom <- expom |> run_sensitivity_analysis(
  base_formula = ~ hs_asthma + hs_child_age + e3_sex,
  methods = c("limma_voom"), scaling_methods = c("none"),
  min_counts_range = c(1,5), min_proportion_range = c(0.1,0.3),
  covariates_to_remove = c("hs_child_age", "e3_sex"),
  pval_col = "adj.P.Val", logfc_col = "logFC",
  pval_threshold = 0.05, stability_metric = "stability_score",
  bootstrap_n = 3, action = "add")
expom |> plot_sensitivity_summary(stability_score_thresh = 0.25, stability_metric = "stability_score")

# Multi-omics integration
expom <- expom |> run_multiomics_integration(method = "MCIA", action = "add")
expom |> plot_factor_summary()

# Factor-outcome association
expom <- expom |> run_association(
  source = "factors",
  outcome = "hs_asthma",
  feature_set = exp_vars[!exp_vars %in% c("hs_asthma", "hs_child_age", "e3_sex")],
  covariates = c("hs_child_age", "e3_sex"),
  action = "add",
  family = "binomial")

# Extract and visualize top contributing features to selected factors
expom <- expom |> extract_top_factor_features(factors = c("V1", "V6", "V9"), method = "percentile", percentile = 0.9, action = "add")
expom |> plot_top_factor_features(top_n = 15, factors = c("V1", "V6", "V9"))

# Find overlapping features across factors
expom <- expom |> run_factor_overlap(stability_score = 0.25, robust_comparison = TRUE, score_col = "stability_score", pval_thresh = 0.05, pval_col = "p.value", logfc_thresh = log2(1.1))
expom |> plot_factor_overlap()

# Correlate DEGs with exposures
expom <- expom |> run_correlation(feature_type = "degs", exposure_cols = exp_vars, action = "add", correlation_cutoff = 0.01, deg_logfc_col = "logFC", deg_logfc_thresh = log2(1), pval_cutoff = 0.1, robust = TRUE, score_col = "stability_score", score_thresh = 0.25)
expom |> plot_correlation_summary(mode = "summary")
expom |> plot_circos_correlation(feature_type = "degs", shared_cutoff = 10, midpoint = 200)

# Network analysis
expom <- expom |> run_create_network(feature_type = "degs", action = "add")
expom |> plot_network(network = "degs", top_n_nodes = 50, include_stats = TRUE, cor_thresh = 0.090, node_color_var = "group", label = TRUE, label_top_n = 5)

# Exposure-omics impact analysis
expom <- expom |> run_exposure_impact(feature_type = "degs")
expom |> plot_exposure_impact(feature_type = "degs", min_per_group = 10, ncol = 2, widths = c(3,1))

# Functional enrichment
expom <- expom |> run_enrichment(geneset = "deg_exp_cor", feature_col = "gene", clustering_approach = "dynamic", pval_threshold = 0.05, pval_col = "p.value", pvalueCutoff = 0.1, pAdjustMethod = "none", qvalueCutoff = 1, logfc_threshold = log2(1), action = "add")
expom |> plot_dotplot_enrichment(geneset = "deg_exp_cor", top_n = 5, n_per_group = 5, add_top_genes = TRUE, top_n_genes = 5)
expom |> plot_go_group_exposures(go_groups = c("Group_1"), feature_col = "gene")

# Pivoting for custom analysis
expom |> pivot_sample() |> head()
expom |> pivot_sample() |> group_by(hs_asthma, e3_sex) |> summarise(n = n())
expom |> pivot_feature() |> head()
expom |> pivot_feature() |> group_by(.exp_name) |> summarise(n = n())
expom |> pivot_exp(omics_name = "Proteomics", features = "TAF7") |> head()
expom |> pivot_exp(omics_name = "Proteomics", features = "TAF7") |>
  ggplot(aes(x = hs_asthma, y = log2(counts), color = hs_asthma, fill = hs_asthma)) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(alpha=0.1) +
  ggpubr::geom_pwc(label = "{p.adj.format}{p.adj.signif}") +
  theme_minimal() +
  ggpubr::rotate_x_text(angle = 45) +
  scale_color_tidy_exp() +
  scale_fill_tidy_exp() +
  labs(x = "", y = expression(Log[2]*"Abd."), fill = "Asthma Status", color = "Asthma Status")

# Final pipeline summary
expom |> run_pipeline_summary(console_print = TRUE, include_notes = TRUE)
