
pop_exposures <- a@metadata$omics_exposure_deg_correlation |> 
  group_by(exposure,assay_name) |> 
  dplyr::reframe(n_assay_name=length(assay_name)) |> 
  group_by(exposure) |> 
  mutate(total=sum(n_assay_name)) |> 
  ungroup() |> 
  filter(exposure %in% c(
    a@metadata$omics_exposure_deg_correlation |> 
      tabyl(exposure) |> 
      arrange(desc(n)) |>
      slice_head(n=30) |> 
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
  theme_pubr(legend="right")+
  theme(plot.title = element_text(face = "bold.italic"))+
  labs(title = "",
       y = "",
       x = "No. of Feature Associations",
       fill = "Assay Name"
  ) 


pop_features <- a@metadata$omics_exposure_deg_correlation |> 
  group_by(feature,category) |> 
  dplyr::reframe(n_category=length(category)) |> 
  group_by(feature) |> 
  mutate(total=sum(n_category)) |> 
  ungroup() |>
  filter(feature %in% c(
    a@metadata$omics_exposure_deg_correlation |> 
      tabyl(feature) |> 
      arrange(desc(n)) |>
      slice_head(n=30) |> 
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
  theme_pubr(legend="right")+
  theme(plot.title = element_text(face = "bold.italic"))+
  labs(title = "",
       y = "",
       x = "No. of Feature Associations",
       fill = "Exposure Category"
  )


.plot_circular_bar(expom@metadata$omics_exposure_deg_correlation |> tabyl(category),"category","n")

.plot_circular_bar(expom@metadata$omics_exposure_deg_correlation |> tabyl(assay_name),"assay_name","n")
                                                                                                                
                               

                                                                                 
                                                                                                    