
# Function for frequently associated exposures
plot_pop_exposures <- function(
    expOmicSet,
    geneset = "degs",
    top_n = 15) {
  require(tidyverse)
  require(janitor)
  require(patchwork)
  require(ggsci)
  require(ggpubr)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  exp_feature_cor_df |> 
      group_by(exposure,assay_name) |>
      dplyr::reframe(n_assay_name=length(assay_name)) |>
      group_by(exposure) |>
      mutate(total=sum(n_assay_name)) |>
      ungroup() |>
      filter(exposure %in% c(
        exp_feature_cor_df |>
          tabyl(exposure) |>
          arrange(desc(n)) |>
          slice_head(n=15) |>
          pull(exposure)
      )) |>
      ggplot(aes(
        x=n_assay_name,
        y=fct_reorder(exposure, total),
        fill=assay_name
      ))+
      geom_bar(stat = "identity",alpha=0.7) +
      scale_fill_npg()+
      scale_color_npg()+
      theme_pubr(legend="right",base_size = 10)+
      theme(plot.title = element_text(face = "bold.italic"))+
      labs(title = "Frequently Associated Exposures",
           y = "",
           x = "No. of Associations",
           fill = "Assay Name"
      )
}

# Function for frequently associated features
plot_pop_features <- function(
    expOmicSet,
    geneset = "degs",
    top_n = 15) {
  require(tidyverse)
  require(janitor)
  require(patchwork)
  require(ggsci)
  require(ggpubr)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  exp_feature_cor_df |> 
      group_by(feature,category) |>
      dplyr::reframe(n_category=length(category)) |>
      group_by(feature) |>
      mutate(total=sum(n_category)) |>
      ungroup() |>
      filter(feature %in% c(
        exp_feature_cor_df |>
          tabyl(feature) |>
          arrange(desc(n)) |>
          slice_head(n=15) |>
          pull(feature)
      )) |>
      ggplot(aes(
        x=n_category,
        y=fct_reorder(feature, total),
        fill=category
      ))+
      geom_bar(stat = "identity",alpha=0.7) +
      scale_fill_cosmic()+
      scale_color_cosmic()+
      theme_pubr(legend="right",base_size = 10)+
      theme(plot.title = element_text(face = "bold.italic"))+
      labs(title = "Frequently Associated Features",
           y = "",
           x = "No. of Associations",
           fill = "Exposure Category"
      )
}

# Function for exposure category associations
plot_exp_assoc <- function(
    expOmicSet,
    geneset = "degs") {
  require(tidyverse)
  require(janitor)
  require(patchwork)
  require(ggsci)
  require(ggpubr)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  exp_feature_cor_df |> 
      tabyl(category) |>
      ggplot(aes(
        x=n,
        y=reorder(category,n),
        fill=category
      ))+
      geom_bar(stat = "identity",alpha=0.7) +
      geom_segment(aes(
        x = n,
        xend = n,
        y = as.numeric(reorder(category, n)) - 0.45,
        yend = as.numeric(reorder(category, n)) + 0.45,
        color = category,
      ), size = 1) +
      scale_fill_cosmic()+
      scale_color_cosmic()+
      theme_pubr(legend="none",base_size = 10)+
      theme(plot.title = element_text(face = "bold.italic"))+
      labs(
           y = "",
           x = "No. of Associations",
           fill = "Exposure Category",
           title = "Exposure Associations")
}

# Function for omic associations
plot_omic_assoc <- function(
    expOmicSet,
    geneset = "degs") {
  require(tidyverse)
  require(janitor)
  require(patchwork)
  require(ggsci)
  require(ggpubr)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  exp_feature_cor_df |> 
      tabyl(assay_name) |>
      ggplot(aes(
        x=n,
        y=reorder(assay_name,n),
        fill=assay_name
      ))+
      geom_bar(stat = "identity",alpha=0.7) +
      geom_segment(aes(
        x = n,
        xend = n,
        y = as.numeric(reorder(assay_name, n)) - 0.45,
        yend = as.numeric(reorder(assay_name, n)) + 0.45,
        color = assay_name,
      ), size = 1) +
      scale_fill_npg()+
      scale_color_npg()+
      theme_pubr(legend="none",base_size = 10)+
      theme(plot.title = element_text(face = "bold.italic"))+
      labs(
           y = "",
           x = "No. of Associations",
           fill = "Exposure Category",
           title = "Omics Associations")
}


plot_exp_feature_assoc_summary <- function(
    expOmicSet,
    geneset = "degs",
    top_n = 15) {
  require(tidyverse)
  require(janitor)
  require(patchwork)
  require(ggsci)
  require(ggpubr)
  
  if(geneset=="degs"){
    if(!"omics_exposure_deg_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_degs()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_deg_correlation
    
  }else if(geneset=="factors"){
    if(!"omics_exposure_factor_correlation" %in% names(expOmicSet@metadata)){
      stop("Please run `correlate_exposures_with_factors()` first.")
    }
    exp_feature_cor_df <- expOmicSet@metadata$omics_exposure_factor_correlation
  } else{
    stop("`geneset` must be either 'degs' or 'factors'")
  }
  
  pop_exposures <- plot_pop_exposures(expOmicSet,
                                      geneset = geneset,
                                      top_n)
  pop_features <- plot_pop_features(expOmicSet,
                                    geneset = geneset,
                                    top_n)
  exp_assoc <- plot_exp_assoc(expOmicSet,
                              geneset = geneset)
  omic_assoc <- plot_omic_assoc(expOmicSet,
                                geneset = geneset)
  
  # Combine using patchwork
  (omic_assoc / exp_assoc) | (pop_features) | (pop_exposures)
}

# require(janitor)
# pop_exposures <- expOmicSet@metadata$omics_exposure_deg_correlation |>
#   group_by(exposure,assay_name) |>
#   dplyr::reframe(n_assay_name=length(assay_name)) |>
#   group_by(exposure) |>
#   mutate(total=sum(n_assay_name)) |>
#   ungroup() |>
#   filter(exposure %in% c(
#     expOmicSet@metadata$omics_exposure_deg_correlation |>
#       tabyl(exposure) |>
#       arrange(desc(n)) |>
#       slice_head(n=15) |>
#       pull(exposure)
#   )) |>
#   ggplot(aes(
#     x=n_assay_name,
#     y=fct_reorder(exposure, total),
#     fill=assay_name
#   ))+
#   geom_bar(stat = "identity",alpha=0.7) +
#   scale_fill_npg()+
#   scale_color_npg()+
#   theme_pubr(legend="right",base_size = 10)+
#   theme(plot.title = element_text(face = "bold.italic"))+
#   labs(title = "Frequently Associated Exposures",
#        y = "",
#        x = "No. of Associations",
#        fill = "Assay Name"
#   )
# pop_exposures
# 
# pop_features <- expOmicSet@metadata$omics_exposure_deg_correlation |>
#   group_by(feature,category) |>
#   dplyr::reframe(n_category=length(category)) |>
#   group_by(feature) |>
#   mutate(total=sum(n_category)) |>
#   ungroup() |>
#   filter(feature %in% c(
#     expOmicSet@metadata$omics_exposure_deg_correlation |>
#       tabyl(feature) |>
#       arrange(desc(n)) |>
#       slice_head(n=15) |>
#       pull(feature)
#   )) |>
#   ggplot(aes(
#     x=n_category,
#     y=fct_reorder(feature, total),
#     fill=category
#   ))+
#   geom_bar(stat = "identity",alpha=0.7) +
#   scale_fill_cosmic()+
#   scale_color_cosmic()+
#   theme_pubr(legend="right",base_size = 10)+
#   theme(plot.title = element_text(face = "bold.italic"))+
#   labs(title = "Frequently Associated Features",
#        y = "",
#        x = "No. of Feature Associations",
#        fill = "Exposure Category"
#   )
# pop_features
# 
# exp_assoc <- expOmicSet@metadata$omics_exposure_deg_correlation |> 
#   tabyl(category) |> 
#   ggplot(aes(
#     x=n,
#     y=reorder(category,n),
#     fill=category
#   ))+
#   geom_bar(stat = "identity",alpha=0.7) +
#   geom_segment(aes(
#     x = n,                    
#     xend = n,                    
#     y = as.numeric(reorder(category, n)) - 0.45,
#     yend = as.numeric(reorder(category, n)) + 0.45,
#     color = category,
#   ), size = 1) +
#   scale_fill_cosmic()+
#   scale_color_cosmic()+
#   theme_pubr(legend="none",base_size = 10)+
#   theme(plot.title = element_text(face = "bold.italic"))+
#   labs(
#        y = "",
#        x = "No. of Associations",
#        fill = "Exposure Category",
#        title = "Exposure Associations")
# 
# exp_assoc
# 
# omic_assoc <- expOmicSet@metadata$omics_exposure_deg_correlation |> 
#   tabyl(assay_name) |> 
#   ggplot(aes(
#     x=n,
#     y=reorder(assay_name,n),
#     fill=assay_name
#   ))+
#   geom_bar(stat = "identity",alpha=0.7) +
#   geom_segment(aes(
#     x = n,                    
#     xend = n,                    
#     y = as.numeric(reorder(assay_name, n)) - 0.45,
#     yend = as.numeric(reorder(assay_name, n)) + 0.45,
#     color = assay_name,
#   ), size = 1) +
#   scale_fill_npg()+
#   scale_color_npg()+
#   theme_pubr(legend="none",base_size = 10)+
#   theme(plot.title = element_text(face = "bold.italic"))+
#   labs(
#        y = "",
#        x = "No. of Associations",
#        fill = "Exposure Category",
#        title = "Omics Associations")
# omic_assoc
# 
# (omic_assoc/ exp_assoc)|pop_features|pop_exposures
# 
# .plot_circular_bar(expOmicSet@metadata$omics_exposure_deg_correlation |> tabyl(category),"category","n")
# 
# .plot_circular_bar(expOmicSet@metadata$omics_exposure_deg_correlation |> tabyl(assay_name),"assay_name","n")




