# --- Plot Functional Enrichment DotPlot: Exp:Feature Correlation ------
plot_exp_enrich_dotplot <- function(
    expOmicSet,
    geneset,
    top_n=10,
    n_per_group=5,
    clustering_approach = "diana"
){
  
  if(!"functional_enrichment" %in% names(expOmicSet@metadata)){
    stop("Please run `run_functional_enrichment` first.")
  }
  
  message("Determining Number of GO Term Clusters...")
  
  go_groups <- expOmicSet@metadata$functional_enrichment[[geneset]] |>
    (\(df) split(df,df$Description) )() |>
    map(~.x |>
          pull(geneID) |>
          paste(collapse ="/") |>
          str_split("/") |>
          unlist() |>
          str_trim() |>
          unique()) |>
    .get_pairwise_overlaps() |>
    dplyr::select(source,target,jaccard) |>
    pivot_wider(names_from = "target",
                values_from = "jaccard") |>
    column_to_rownames("source") |>
    as.matrix() |>
    .cluster_mat(clustering_approach = clustering_approach) |>
    (\(x) df=data.frame(
      Description=names(x),
      go_group=as.numeric(x)))()
  
  go_group_df <- expOmicSet@metadata$functional_enrichment[[geneset]] |>
    left_join(go_groups,
              by="Description")
  
  message("Plotting Enrichment Results...")
  go_group_df |>
    filter(go_group %in% c(
      go_group_df |> 
        group_by(go_group) |>
        reframe(
          score=mean(-log10(p.adjust) * Count)) |> 
        arrange(desc(score)) |> 
        slice_head(n=top_n) |> 
        pull(go_group)
    )) |> 
    group_by(category,exp_name,go_group) |>
    arrange(p.adjust) |>
    slice_head(n=n_per_group) |>
    ungroup() |>
    ggplot(aes(
      x = fct_reorder(exp_name, -Count), 
      y = fct_reorder(Description, -log10(p.adjust)),  
      size = Count,
      color = -log10(p.adjust)
    )) +
    geom_point() +
    theme_pubr(legend = "right",
               base_size = 10) +
    facet_grid(go_group~category, space ="free",scales = "free") +
    rotate_x_text(angle=45)+
    scale_color_gradient(low="thistle1",high="magenta4") +
    theme(strip.text = element_text(face = "bold.italic",angle=90))+
    labs(
      x="",
      y="",
      title="",
      size="Number of Genes",
      color=expression(-Log[10](P)),
    )
}




a@metadata$functional_enrichment$deg_exp_cor |>
  #filter(category == "indoor_pollution") |>
  group_by(exp_name,category) |>
  arrange(p.adjust) |>
  slice_head(n=5) |>
  ungroup() |>
  ggplot(aes(
    x=exp_name,
    y=fct_reorder(Description,as.numeric(as.factor(Description))),
    color=-log10(p.adjust),
    size=Count
  )) +
  geom_point() +
  theme_pubr(legend = "right",
             base_size = 10) +
  facet_grid(direction~category, scales="free_y") +
  rotate_x_text(angle=45)+
  scale_color_gradient(low="thistle1",high="magenta4") +
  theme(strip.text = element_text(face = "bold.italic"))+
  labs(
    x="",
    y="",
    title="",
    size="Number of Genes",
    color=expression(-Log[10](P)),
  )


go_groups <- a@metadata$functional_enrichment$deg_exp_cor |>
  (\(df) split(df,df$Description) )() |>
  map(~.x |>
        pull(geneID) |>
        paste(collapse ="/") |>
        str_split("/") |>
        unlist() |>
        str_trim() |>
        unique()) |>
  .get_pairwise_overlaps() |>
  dplyr::select(source,target,jaccard) |>
  pivot_wider(names_from = "target",
              values_from = "jaccard") |>
  column_to_rownames("source") |>
  as.matrix() |>
  .cluster_mat(clustering_approach = "diana") |>
  (\(x) df=data.frame(Description=names(x),cluster=as.numeric(x)))()

d <- a@metadata$functional_enrichment$deg_exp_cor |>
  left_join(go_groups,by="Description")

a@metadata$functional_enrichment$deg_exp_cor |>
  left_join(go_groups,by="Description") |>
  mutate(cluster=factor(cluster)) |>
  arrange(cluster) |>
  #mutate(cluster=fct_reorder(cluster,cluster)) |>
  ggplot(aes(
    x=exp_name,
    y=fct_reorder(Description,as.numeric(as.factor(Description))),
    color=cluster,
    size=Count
  )) +
  geom_point() +
  theme_pubr(legend = "right",
             base_size = 10) +
  facet_grid(.~category, scales="free_y") +
  rotate_x_text(angle=45)+
  scale_color_brewer(palette="Set1") +
  theme(strip.text = element_text(face = "bold.italic"))+
  labs(
    x="",
    y="",
    title="",
    size="Number of Genes",
    color="Cluster",
  )

d |>
  group_by(category,exp_name,cluster) |>
  arrange(p.adjust) |>
  slice_head(n=1) |>
  ungroup() |>
  group_by(Description) |>
  mutate(max_n=max(Count)) |>
  ungroup() |>
  arrange(desc(max_n)) |>
  mutate(Description=fct_reorder(Description,max_n)) |>
  ggplot(aes(
    x=Count,
    y=Description,
    color=category
  ))+
  geom_point()+
  geom_segment(aes(
    x=0,
    xend=Count,
    y=Description,
    yend=Description
  ))+
  #facet_grid(cluster~category,space = "free",scales = "free")+
  theme_pubr()


d |>
  filter(cluster %in% c(
    d |> 
      group_by(cluster) |>
      reframe(score=mean(-log10(p.adjust)*Count)) |> 
      arrange(desc(score)) |> 
      slice_head(n=5) |> 
      pull(cluster)
  )) |> 
  group_by(category,exp_name,cluster) |>
  arrange(p.adjust) |>
  slice_head(n=5) |>
  ungroup() |>
  ggplot(aes(
    x = fct_reorder(exp_name, -Count),  # Order assays by number of genes
    y = fct_reorder(Description, -log10(p.adjust)),  # Order GO terms by significance
    size = Count,
    color = -log10(p.adjust)
  )) +
  geom_point() +
  theme_pubr(legend = "right",
             base_size = 10) +
  facet_grid(cluster~category, space ="free",scales = "free") +
  rotate_x_text(angle=45)+
  scale_color_gradient(low="thistle1",high="magenta4") +
  theme(strip.text = element_text(face = "bold.italic",angle=90))+
  labs(
    x="",
    y="",
    title="",
    size="Number of Genes",
    color=expression(-Log[10](P)),
  )

# --- Gene Plot Per Cluster ------------

x=d |>
  filter(cluster %in% c(
    d |> 
      group_by(cluster) |>
      reframe(score=mean(-log10(p.adjust)*Count)) |> 
      arrange(desc(score)) |> 
      slice_head(n=5) |> 
      pull(cluster)
  )) |> 
  group_by(category,exp_name,cluster) |>
  arrange(p.adjust) |>
  slice_head(n=5) |>
  ungroup() |> 
  separate_rows(geneID,sep="/") |> 
  dplyr::select(exp_name,Description,geneID,cluster) |> 
  dplyr::group_split(cluster) |>
  map( ~ .x |> 
         tabyl(geneID))
  
# --- Test ---------

d |>
  group_by(category,exp_name,cluster) |>
  arrange(p.adjust) |>
  slice_head(n=5) |>
  ungroup() |>
  filter(cluster==2) |>
  ggplot(aes(
    x = fct_reorder(exp_name, -Count),  # Order assays by number of genes
    y = fct_reorder(Description, -log10(p.adjust)),  # Order GO terms by significance
    #size = num_genes,
    fill = -log10(p.adjust)
  )) +
  geom_tile() +
  theme_pubr(legend = "right",
             base_size = 10) +
  facet_grid(cluster~category, space ="free",scales = "free") +
  rotate_x_text(angle=45)+
  scale_fill_gradient(low="thistle1",high="magenta4") +
  theme(strip.text = element_text(face = "bold.italic",angle = 90))+
  labs(
    x="",
    y="",
    title="",
    size="Number of Genes",
    fill=expression(-Log[10](P)),
  )




e <- d |>
  (\(df) split(df,df$cluster))() |>
  map(~.x |>
        pull(genes) |>
        paste(collapse = ", ") |>
        str_split(",") |>
        unlist() |>
        str_trim() |>
        unique()) |>
  map(~ expom@metadata$sensitivity_analysis$sensitivity_df |>
        filter(molecular_feature %in% .x) |>
        filter(adj.P.Val<0.05) |>
        group_by(molecular_feature) |>
        dplyr::summarize(
          mean_logfc=mean(abs(logFC))
        )) |>
  bind_rows(.id="cluster")


e |>
  group_by(cluster) |>
  dplyr::mutate(
    mean_group=mean(mean_logfc)
  ) |>
  ggplot(aes(
    x=reorder(cluster,-mean_group),
    y=mean_logfc
  ))+
  geom_violin()+
  theme_pubr()+
  labs(
    x="Cluster",
    y="Mean Abs. LogFC",
    title=""
  )
