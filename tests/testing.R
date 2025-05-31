


a=expom_multi

var_map <- a@metadata$var_info |>
  dplyr::select(variable) |>
  mutate(exp_name="exposure") |>
  bind_rows(a@metadata$common_top_factor_features |>
              ungroup() |>
              dplyr::select(feature,exp_name) |>
              setNames(c("variable","exp_name")))

min_map <- var_map |>
  sample_n(100) |>
  as.data.frame()

a=a |>
  correlate_exposoures_omics(
    variable_map = min_map,
    correlation_method = "spearman",
    correlation_cutoff = 0.03,
    cor_pval_column = "p.value",
    pval_cutoff = 0.05,
    action = "add"
  )


b=MultiAssayExperiment::colData(a) |>
  as.data.frame() |>
  dplyr::select(
    all_of(c(exp_vars,"hs_asthma"))) |>
  select_if(is.numeric) |>
  mutate(hs_asthma = factor(hs_asthma,
                            levels = c("non_asthmatic", "asthma"))

model <- glm(hs_asthma ~ ., data = b, family = "binomial")
