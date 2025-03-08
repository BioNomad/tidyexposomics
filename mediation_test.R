library(dplyr)
library(purrr)
library(stringr)

# Extract enrichment results and merge with Go Groups
merged_enrich_res <- MultiAssayExperiment::metadata(expom_enrich)$functional_enrichment |> 
  pluck("deg_exp_cor") |> 
  pluck("enrich_res") |> 
  inner_join(
    MultiAssayExperiment::metadata(expom_enrich)$functional_enrichment |> 
      pluck("deg_exp_cor") |> 
      pluck("go_groups"), 
    by = "Description") |> 
  mutate(go_group=paste("Group",go_group,sep="_"))

# --- New Test ------
outcome <- "pftfev1fvc_actual"
exposure <- "pm25f"
covariates <- c("age","sex","race")



a <- merged_enrich_res |> 
  dplyr::select(exp_name,Cluster,go_group) |> 
  distinct() |> 
  filter(!exp_name %in% c("CD4+ T-cell miRNA","CD16+ Monocyte miRNA")) |> 
  slice_head(n=5)

b <- a |>
  pmap(function(exp_name, Cluster, go_group) {
    df_filt <- merged_enrich_res |> 
      filter(exp_name == !!exp_name,
             Cluster == !!Cluster, 
             go_group == !!go_group)
    
    genes_per_go_group <- df_filt |> 
      pull(geneID) |> 
      str_split("/") |> 
      unlist() |> 
      unique()
    
    df <- df_filt |> 
      dplyr::select(exp_name, Cluster, go_group) |>
      distinct() |> 
      mutate(genes=paste(genes_per_go_group, collapse=",")) 
    
    exp <- .update_assay_colData(expom_enrich,exp_name)
    
    assay <- assay(exp)
    
    assay_filt <- assay[rownames(assay) %in% genes_per_go_group,] |> 
      t()
    
    assay_filt <- log2(assay_filt + 1)
    
    pca_res <- data.frame(PC_exp_mat=prcomp(assay_filt,scale. = TRUE)$x[,1],
                      id_to_map=names(prcomp(assay_filt,scale. = TRUE)$x[,1]))
    
    merged <- colData(exp) |> 
      as.data.frame() |> 
      rownames_to_column("id_to_map") |>
      inner_join(pca_res, 
                 by="id_to_map") 
    
    merged <- merged |> 
      dplyr::select(all_of(c("PC_exp_mat",outcome,exposure,covariates))) |> 
      drop_na()
      
    # Run Mediation Analysis
    model_total <- lm(
      as.formula(paste(outcome,"~",exposure,"+",
                       paste(covariates,collapse = "+"))), 
      data = merged)
    
    model_mediator <- lm(
      as.formula(paste("PC_exp_mat ~",exposure,"+",
                       paste(covariates,collapse = "+"))), 
      data = merged)
    
    model_outcome <- lm(
      as.formula(paste(outcome, "~ PC_exp_mat + ",exposure,"+",
                       paste(covariates,collapse = "+"))), 
      data = merged)
    
    med_res <- mediate(
      model_mediator, 
      model_outcome, 
      covariates = covariates,
      treat = "pm25f", 
      mediator = "PC_exp_mat", 
      boot = TRUE, 
      sims = 10)
    
    return(med_res)
  }) 

# --- Med Test -------
outcome <- "fev_height"
exposure <- "pm25f"
covariates <- c("age","sex","race","bmipercentile","income5","smoking","cbceosinnum","pm10f","no2f","sige_total_ige","urine_Cu")
others <- c("pftprefev1best","gli_height","pftfev1fvc_actual","pm10f")


merged <- aw_dataset_cv2_filt |> 
  dplyr::select(all_of(c(exposure,covariates,others))) |> 
  mutate(smoking=ifelse(smoking==1,"smoking","non_smoking")) |>
  drop_na() |> 
  mutate_if(is.numeric, ~ as.numeric(scale(.x)))  

merged=merged |> mutate_if(is.factor,~as.character(.))

# Run Mediation Analysis
model_total <- lm(
  as.formula(paste(outcome,"~",exposure,"+",
                   paste(covariates,collapse = "+"))), 
  data = merged)


merged |> glimpse()

model_mediator <- lm(
  as.formula(paste("PC_exp_mat ~",exposure,"+",
                   paste(covariates,collapse = "+"))), 
  data = merged)

model_outcome <- lm(
  as.formula(paste(outcome, "~",exposure,"+", "PC_exp_mat +",
                   paste(covariates,collapse = "+"))), 
  data = merged)


model_mediator <- glm(
  as.formula(paste("PC_exp_mat ~",exposure)), 
  data = merged)

model_outcome <- glm(
  as.formula(paste(outcome, "~", exposure,"+","PC_exp_mat"
                   )), 
  data = merged)

med_res <- mediate(
  model_mediator, 
  model_outcome, 
  treat = "pm25f", 
  mediator = "PC_exp_mat", 
  boot = TRUE, 
  sims = 100)
# --- New New Test ---------
merged$race <- factor(merged$race)  # Might cause issues due to extreme imbalance
merged$smoking <- factor(merged$smoking)
merged$sex <- factor(merged$sex)
merged$income5 <- factor(merged$income5)
merged <- merged |> 
  drop_na()

outcome <- "fev_height"
exposure <- "pm25f"
covariates <- c("age","sex")
others <- c("pftprefev1best","gli_height","pftfev1fvc_actual","pm10f","cbceosinnum","race")


merged <- aw_dataset_cv2_filt |> 
  dplyr::select(all_of(c(exposure,covariates,others,outcome))) |> 
  #mutate(smoking=ifelse(smoking==1,"smoking","non_smoking")) |>
  drop_na() |> 
  mutate_if(is.numeric, ~ as.numeric(scale(.x)))  |> 
  mutate_if(is.character,~ as.factor(.))

# Re-run models with corrected covariates
mediator_formula <- as.formula(
  paste("cbceosinnum ~ pm25f +", paste(covariates, collapse = " + "))
)

med_res <- mediate(mediation_results[[1]],
                   mediation_results[[2]],
                   treat = "exposure_var",
                   mediator = "mediator_var",
                   boot = TRUE,
                   sims = 1000)

outcome_formula <- as.formula(
  paste("fev_height ~ pm25f + cbceosinnum +", paste(covariates, collapse = " + "))
)

# Fit models
mediator_model <- lm(mediator_formula, data = merged)
outcome_model <- lm(outcome_formula, data = merged)

# Mediation analysis
med_res <- mediate(
  model.m = mediator_model, 
  model.y = outcome_model, 
  treat = "pm25f", 
  mediator = "cbceosinnum", 
  boot = TRUE, 
  sims = 1000
)

mediation_results <- .run_mediation(
  data = merged,
  outcome = "fev_height",
  exposure = "pm25f",
  mediator = "cbceosinnum",
  covariates = c("age", "sex", "race")
)
# --- Old Test ------


# Get all unique (exp_name, Cluster) combinations
unique_combinations <- merged_enrich_res |> 
  distinct(exp_name, Cluster)

# Iterate over combinations
a <- unique_combinations |> 
  pmap(function(exp_name, Cluster) {
    df_filt <- merged_enrich_res |> 
      filter(exp_name == !!exp_name, Cluster == !!Cluster) 
    
    genes_per_go_group <- unique(df_filt$go_group) |> 
      set_names() |> 
      map(~ merged_enrich_res |> 
            filter(go_group == .x) |> 
            pull(geneID) |> 
            str_split("/") |> 
            unlist() |> 
            unique()
      )
    
    return(genes_per_go_group)
  })

a <- unique(merged_enrich_res$exp_name) |> 
  set_names() |> 
  map( ~ {
    # Grab results from one omic
    df_exp <- merged_enrich_res |> 
      filter(exp_name == .x)
    
    # Get df per Cluster
    res_per_group_per_cluster <- unique(df_exp$Cluster) |> 
      set_names() |> 
      map( ~ {
        df_go_clust <- df_exp |> 
          filter(Cluster == .x) 
        
        # Get genes per omic, cluster and go group
        res <- unique(df_go_clust$go_group) |> 
          set_names() |> 
          map(~ df_go_clust |> 
                filter(go_group == .x) |> 
                pull(geneID) |> 
                str_split("/") |> 
                unlist() |> 
                unique()
          )
        
        return(res)
      })
    
    return(res_per_group_per_cluster)
    
  })

b <- names(a) |> 
  set_names() |> 
  map(~ {
    exp <- .update_assay_colData(expom_enrich,.x)
    
  })

# For each omic, and then go groups in that omic, create a matrix of feature by 
# samples, get the first principal component
b <- a |> 
  map(~ map(.x, ~ expom_enrich[[.x]])) |> 
  map(~ map(.x, ~ t(.x))) |> 
  map(~ map(.x, ~ prcomp(.x)$x[,1])) |> 
  map(~ do.call(cbind, .x))

# --- Testing GO Group PCA ----

a <- MultiAssayExperiment::metadata(expom_enrich)$functional_enrichment |> 
  pluck("deg_exp_cor") |> 
  pluck("enrich_res") |> 
  inner_join(
    MultiAssayExperiment::metadata(expom_enrich)$functional_enrichment |> 
      pluck("deg_exp_cor") |> 
      pluck("go_groups"), 
    by = "Description") |> 
  mutate(go_group=paste("Group",go_group,sep="_")) |> 
  filter(!grepl("miRNA",exp_name))

b <- a |>
  #filter(go_group %in% c("Group_2")) |> 
  dplyr::select(exp_name, Cluster, go_group) |> 
  distinct() |> 
  pmap(function(exp_name, Cluster, go_group) {
    
    # Filter enrichment results
    df_filt <- a |> 
      filter(exp_name == !!exp_name,
             Cluster == !!Cluster, 
             go_group == !!go_group)
    
    # Get genes per GO group
    genes_per_go_group <- df_filt |> 
      pull(geneID) |> 
      str_split("/") |> 
      unlist() |> 
      unique()
    
    # Get top 3 GO term descriptions
    top_terms <- df_filt |> 
      dplyr::pull(Description) |> 
      unique() |> 
      head(3) |> 
      paste(collapse = "; ")
    
    # If there are less than 10 genes, skip
    if(length(genes_per_go_group) < 10) {
      message(paste("Skipping", exp_name, Cluster, go_group, "- not enough genes"))
      return(NULL)
    }
    
    # Update assay and colData
    exp <- .update_assay_colData(expom_enrich, exp_name)
    
    # Get assay data
    assay <- assay(exp)
    
    # Filter assay for selected genes
    assay_filt <- assay[rownames(assay) %in% genes_per_go_group, , drop = FALSE]
    
    # Ensure assay_filt has at least 2 features and variance
    if (nrow(assay_filt) < 2 || all(apply(assay_filt, 1, var) == 0)) {
      message(paste("Skipping", exp_name, Cluster, go_group, "- insufficient variance in features"))
      return(NULL)
    }
    
    # Transpose to make samples rows
    assay_filt <- t(assay_filt)
    
    # Apply log transformation to stabilize variance
    assay_filt <- log2(assay_filt + 1)
    
    # Perform PCA safely
    pca_res <- prcomp(assay_filt, scale. = TRUE)
    
    # Create output dataframe
    pca_res_df <- data.frame(
      PC_exp_mat = pca_res$x[,1],
      id_to_map = rownames(pca_res$x)
    )
    
    # Rename PC1 column to indicate the source
    names(pca_res_df) <- c(paste("PC", exp_name, Cluster, go_group, sep="/"), "id_to_map")
    
    # Create dataframe for genes and terms
    pca_genes_terms_df <- data.frame(
      term = paste("PC", exp_name, Cluster, go_group, sep="/"),
      top_terms = top_terms,
      genes = paste(genes_per_go_group, collapse = ","),
    )
    
    return(list(pca_res_df=pca_res_df,
                pca_genes_terms_df=pca_genes_terms_df))
  })

c=b[unlist(purrr:::map(b,~{!is.null(.x)}))]

# d=purrr::reduce(c,inner_join, by = "id_to_map") |> 
#   inner_join(colData(expom_enrich) |>
#                as.data.frame() |> 
#                dplyr::select(fev_height) |> 
#                rownames_to_column("id_to_map"),
#              by = "id_to_map") |> 
#   column_to_rownames("id_to_map")

# # correlate data
# cor_res <- cor(d, use = "pairwise.complete.obs") |> 
#   as.data.frame() |> 
#   rownames_to_column("var1") |> 
#   pivot_longer(-var1, names_to = "var2", values_to = "correlation") |> 
#   filter(!is.na(correlation),
#          var1 != var2) |> 
#   #mutate(correlation=abs(correlation)) |> 
#   arrange(correlation) |>
#   filter(var1=="fev_height")
  
# filter list to just the pca_res_df
c_pca_res_df <- purrr::map(c,~.x$pca_res_df)

# filter list to just the pca_genes_terms_df
c_pca_genes_terms_df <- purrr::map(c,~.x$pca_genes_terms_df) |> 
  bind_rows()


e <- purrr::reduce(c_pca_res_df,inner_join, by = "id_to_map") |> 
  inner_join(colData(expom_enrich) |>
               as.data.frame() |> 
               dplyr::select(fev_height,age,sex,race,smoking) |> 
               rownames_to_column("id_to_map"),
             by = "id_to_map") |> 
  column_to_rownames("id_to_map")

#colnames(e) <- gsub(" |\\+|\\.|\\-","_",colnames(e))

pc_cols <- colnames(e)[grepl("PC/",colnames(e))]

e <- e |> 
  mutate(across(c(pc_cols,fev_height), ~ as.numeric(scale(.x))))

pc_assoc_res <- purrr::map(
  pc_cols,
  function(pc_col) {
    # Create a placeholder column
    e$placeholder_col <- e[[pc_col]]
      
    model <- glm(
      as.formula(paste("fev_height ~ placeholder_col+age+sex+race+smoking")),data = e)
    
    return(broom::tidy(model) |> 
             mutate(term=ifelse(term=="placeholder_col",pc_col,term)))
  }
) |> 
  bind_rows() |> 
  filter(grepl("PC/",term)) |> 
  inner_join(c_pca_genes_terms_df, by = c("term"="term")) 


# gam_assoc_res <- purrr::map(
#   pc_cols,
#   function(pc_col) {
#     
#     # Create a placeholder column
#     e$placeholder_col <- e[[pc_col]]
#     
#     model <- gam(
#       as.formula(paste("fev_height ~ s(placeholder_col) + age + sex + race")), # s() applies a smooth function
#       data = e,
#       method = "REML"
#     )
#     return(broom::tidy(model) |> 
#              mutate(term=ifelse(term=="s(placeholder_col)",pc_col,term)))
#   }
# ) |> 
#   bind_rows()
