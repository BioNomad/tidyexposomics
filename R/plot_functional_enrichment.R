# expom@metadata$exposure_category_enrichment_summary |> 
#   #filter(category == "indoor_pollution") |> 
#   group_by(assay_name,category) |> 
#   arrange(min_p_adj) |> 
#   slice_head(n=5) |> 
#   ungroup() |> 
#   ggplot(aes(
#     x=assay_name,
#     y=fct_reorder(Description,as.numeric(as.factor(Description))),
#     color=-log10(min_p_adj),
#     size=num_genes
#   )) +
#   geom_point() +
#   theme_pubr(legend = "right",
#              base_size = 10) +
#   facet_grid(.~category, scales="free_y") +
#   rotate_x_text(angle=45)+
#   scale_color_gradient(low="thistle1",high="magenta4") +
#   theme(strip.text = element_text(face = "bold.italic"))+
#   labs(
#     x="",
#     y="",
#     title="",
#     size="Number of Genes",
#     color=expression(-Log[10](P)),
#   ) 
# 
# 
# go_groups <- expom@metadata$exposure_category_enrichment_summary |> 
#   (\(df) split(df,df$Description) )() |> 
#   map(~.x |> 
#         pull(genes) |> 
#         paste(collapse =", ") |>
#         str_split(",") |> 
#         unlist() |> 
#         str_trim() |> 
#         unique()) |> 
#   .get_pairwise_overlaps() |> 
#   dplyr::select(source,target,jaccard) |> 
#   pivot_wider(names_from = "target",
#               values_from = "jaccard") |> 
#   column_to_rownames("source") |> 
#   as.matrix() |> 
#   .cluster_mat(clustering_approach = "diana") |> 
#   (\(x) df=data.frame(Description=names(x),cluster=as.numeric(x)))()
# 
# 
# go_with_cluster <- expom@metadata$exposure_category_enrichment_summary |>
#   left_join(c,by="Description") |> 
#   mutate(cluster=factor(cluster)) |> 
#   arrange(cluster) |> 
#   mutate(cluster=fct_reorder(cluster,cluster)) |> 
#   ggplot(aes(
#     x=assay_name,
#     y=fct_reorder(Description,as.numeric(as.factor(Description))),
#     color=cluster,
#     size=num_genes
#   )) +
#   geom_point() +
#   theme_pubr(legend = "right",
#              base_size = 10) +
#   facet_grid(.~category, scales="free_y") +
#   rotate_x_text(angle=45)+
#   scale_color_brewer(palette="Set1") +
#   theme(strip.text = element_text(face = "bold.italic"))+
#   labs(
#     x="",
#     y="",
#     title="",
#     size="Number of Genes",
#     color="Cluster",
#   )
# 
# d |> 
#   group_by(category,assay_name,cluster) |> 
#   arrange(min_p_adj) |> 
#   slice_head(n=1) |> 
#   ungroup() |>
#   group_by(Description) |> 
#   mutate(max_n=max(num_genes)) |>
#   ungroup() |> 
#   arrange(desc(max_n)) |> 
#   mutate(Description=fct_reorder(Description,max_n)) |>
#   ggplot(aes(
#     x=num_genes,
#     y=Description,
#     color=category
#   ))+
#   geom_point()+
#   geom_segment(aes(
#     x=0,
#     xend=num_genes,
#     y=Description,
#     yend=Description
#   ))+
#   #facet_grid(cluster~category,space = "free",scales = "free")+
#   theme_pubr()
#   
# 
# d |> 
#   group_by(category,assay_name,cluster) |> 
#   arrange(min_p_adj) |> 
#   slice_head(n=1) |> 
#   ungroup() |> 
#   ggplot(aes(
#     x = fct_reorder(assay_name, -num_genes),  # Order assays by number of genes
#     y = fct_reorder(Description, -log10(min_p_adj)),  # Order GO terms by significance
#     size = num_genes,
#     color = -log10(min_p_adj)
#   )) +
#   geom_point() +
#   theme_pubr(legend = "right",
#              base_size = 10) +
#   facet_grid(cluster~category, space ="free",scales = "free") +
#   rotate_x_text(angle=45)+
#   scale_color_gradient(low="thistle1",high="magenta4") +
#   theme(strip.text = element_text(face = "bold.italic"))+
#   labs(
#     x="",
#     y="",
#     title="",
#     size="Number of Genes",
#     color=expression(-Log[10](P)),
#   ) 
# 
# d |> 
#   group_by(category,assay_name,cluster) |> 
#   arrange(min_p_adj) |> 
#   slice_head(n=5) |> 
#   ungroup() |> 
#   filter(cluster==2) |> 
#   ggplot(aes(
#     x = fct_reorder(assay_name, -num_genes),  # Order assays by number of genes
#     y = fct_reorder(Description, -log10(min_p_adj)),  # Order GO terms by significance
#     #size = num_genes,
#     fill = -log10(min_p_adj)
#   )) +
#   geom_tile() +
#   theme_pubr(legend = "right",
#              base_size = 10) +
#   facet_grid(cluster~category, space ="free",scales = "free") +
#   rotate_x_text(angle=45)+
#   scale_fill_gradient(low="thistle1",high="magenta4") +
#   theme(strip.text = element_text(face = "bold.italic",angle = 90))+
#   labs(
#     x="",
#     y="",
#     title="",
#     size="Number of Genes",
#     fill=expression(-Log[10](P)),
#   ) 
# 
# 
# 
# 
# e <- d |> 
#   (\(df) split(df,df$cluster))() |> 
#   map(~.x |> 
#         pull(genes) |> 
#         paste(collapse = ", ") |> 
#         str_split(",") |> 
#         unlist() |> 
#         str_trim() |> 
#         unique()) |> 
#   map(~ expom@metadata$sensitivity_analysis$sensitivity_df |> 
#         filter(molecular_feature %in% .x) |> 
#         filter(adj.P.Val<0.05) |> 
#         group_by(molecular_feature) |> 
#         dplyr::summarize(
#           mean_logfc=mean(abs(logFC))
#         )) |> 
#   bind_rows(.id="cluster") 
# 
# 
# e |> 
#   group_by(cluster) |>
#   dplyr::mutate(
#     mean_group=mean(mean_logfc)
#   ) |> 
#   ggplot(aes(
#     x=reorder(cluster,-mean_group),
#     y=mean_logfc
#   ))+
#   geom_violin()+
#   theme_pubr()+
#   labs(
#     x="Cluster",
#     y="Mean Abs. LogFC",
#     title=""
#   ) 
