x <- a@metadata$functional_enrichment$deg_exp_cor$enrich_res |>
  left_join(a@metadata$functional_enrichment$deg_exp_cor$go_groups,
            by="Description")

y <- x |>
  filter(go_group %in% c(
    x |> 
      group_by(go_group) |>
      reframe(
        score=mean(-log10(p.adjust) * Count)) |> 
      arrange(desc(score)) |> 
      slice_head(n=5) |> 
      pull(go_group)
  )) |> 
  group_by(category,exp_name,go_group) |>
  arrange(p.adjust) |>
  slice_head(n=5) |>
  ungroup() 


yy <- y |> 
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

z <- y |>
  separate_rows(geneID,sep="/") |> 
  group_by(go_group,geneID) |> 
  reframe(n=n()) |> 
  group_by(go_group) |>
  arrange(desc(n)) |> 
  mutate(
    geneID=fct_reorder(geneID,n)
  ) |>
  slice_head(n=10)


zz <- z |> 
  ggplot(aes(
    x=n,
    y=geneID,
    color=n
  ))+
  geom_point()+
  geom_segment(aes(
    x=0,
    xend=n,
    y=geneID,
    yend=geneID
  ))+
  theme_pubr(legend="none",base_size = 10)+
  scale_color_gradient(low = "skyblue",high = "midnightblue")+
  theme(strip.text = element_text(face="bold.italic"))+
  facet_grid(go_group~.,scales="free_y")+
  labs(
    x="Frequency",
    y=""
  )

yy|zz
