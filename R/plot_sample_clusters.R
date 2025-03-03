plot_sample_clusters <- function(
    expomicset,
    cols_of_interest=NULL
){
  require(ggpubr)
  require(patchwork)
  require(broom)
  require(tidyverse)
  require(grid)
  require(gridExtra)
  require(ComplexHeatmap)
  require(MultiAssayExperiment)
  require(SummarizedExperiment)
  
  if(!"sample_clustering" %in% names(expomicset@metadata)){
    stop("Please run `cluster_samples()` first")
  }
  
  if(is.null(cols_of_interest)){
    cols_of_interest <- colnames(expomicset)
  } else{
    cols_of_interest <- cols_of_interest[cols_of_interest %in% colnames(colData(expomicset))]
  }
  
  # add in sample group to the colData
  meta <- expomicset |> 
    colData() |> 
    as.data.frame() |> 
    rownames_to_column("id_to_map") |>
    left_join(
      data.frame(
        id_to_map = names(expomicset@metadata$sample_clustering$sample_groups),
        cluster = as.vector(expomicset@metadata$sample_clustering$sample_groups)
      ),
      by="id_to_map"
    ) |> 
    dplyr::select(all_of(c(cols_of_interest,"cluster"))) |> 
    mutate_at(vars(cols_of_interest), as.numeric)
  
  # perform tests to determine what variables are higher or lower in each cluster using anova if more than 2 clusters and t-test if two clusters
  if(length(unique(meta$cluster)) > 2){
    tests <- lapply(
      cols_of_interest,
      function(col){
        aov_res <- aov(
          as.formula(paste0(col, " ~ cluster")),
          data=meta
        ) |> 
          broom::tidy() |> 
          mutate(variable=col)
      }
    ) |> 
      bind_rows() |> 
      filter(!is.na(p.value))
  } else {
    tests <- lapply(
      cols_of_interest,
      function(col){
        ttest_res <- t.test(
          meta |> filter(cluster==1) |> pull(!!sym(col)),
          meta |> filter(cluster==2) |> pull(!!sym(col))
        ) |> 
          broom::tidy() |> 
          dplyr::mutate(variable=col)
      }
    ) |> 
      bind_rows()
  }
  
  lollipop_plot <- tests |> 
    mutate(sig=case_when(
      p.value < 0.05 & statistic > 0 ~ "up",
      p.value < 0.05 & statistic < 0 ~ "down",
      .default = "ns"
    )) |>
    ggplot(aes(x=statistic,
               y=reorder(variable,statistic),
               color=sig)) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_segment(aes(
      x=0,
      y=variable,
      xend=statistic
    ))+
    geom_point() +
    theme_minimal() +
    scale_color_manual(values = c(
      "up" = "red3",
      "down" = "blue4",
      "ns" = "grey55"
    ))+
    theme(legend.position = "none")+
    labs(
      x="Test-Statistic",
      y="",
      title=paste("Cluster Differences in", "Variables of Interest",sep="\n"),
    )
  
  wrap_plots(lollipop_plot,
             grid.grabExpr(
               draw(expomicset@metadata$sample_clustering$heatmap)))+
    plot_layout(widths = c(1,3))
  
}





