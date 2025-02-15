plot_normality_sum <- function(
    expOmicSet
){
  
  norm_df <- expOmicSet@metadata$normality$norm_df
  
  # Create normality plot
  norm_plot <- table(norm_df$p.value > 0.05) |>
    as.data.frame() |>
    mutate(Var1 = case_when(
      Var1 == "FALSE" ~ "Not Normal",
      Var1 == "TRUE" ~ "Normal"
    )) |>
    ggplot(aes(
      x = Var1,
      y = Freq,
      fill = Var1
    )) +
    geom_bar(stat = "identity", alpha = 0.5) +
    geom_segment(
      aes(x = as.numeric(as.factor(Var1)) - 0.45,
          xend = as.numeric(as.factor(Var1)) + 0.45,
          y = Freq,
          yend = Freq,
          color = Var1),
      size = 1
    ) +
    theme_pubr(legend = "right") +
    scale_fill_lancet() +
    scale_color_lancet(guide = FALSE) +
    labs(
      x = "",
      y = "No. of Exposures",
      fill = "",
      title = "Normality of Exposure Variables",
      subtitle = "Shapiro-Wilk Test"
    )
  
  transform_plot <- expOmicSet@metadata$transformation$norm_summary |> 
    mutate(max_normal=ifelse(normal==max(normal),1,0)) |> 
    dplyr::select(transformation, normal,not_normal,max_normal) |>
    pivot_longer(cols = -c(transformation,max_normal),
                 names_to = "normal",
                 values_to = "Freq") |>
    arrange(desc(max_normal)) |> 
    mutate(transformation=case_when(
      transformation=="log_trans" ~ "Log2",
      transformation=="sqrt_trans" ~ "Square Root",
      transformation=="x_1_3_trans" ~ "X^(1/3)",
      #transformation=="x_1_3_trans" ~ expression(X^{1/3}),
      transformation=="none_trans" ~ "No Transformation"
    )) |>
    mutate(normal=case_when(
      normal=="normal" ~ "Normal",
      normal=="not_normal" ~ "Not Normal"
    )) |> 
    mutate(transformation=factor(transformation)) |> 
    ggplot(aes(
      x = reorder(transformation,-max_normal),
      y = Freq,
      fill = normal
    )) +
    geom_bar(stat = "identity", alpha = 0.5) +
    theme_pubr(legend = "right") +
    rotate_x_text(angle=45)+
    scale_fill_lancet() +
    scale_color_lancet(guide = FALSE) +
    labs(
      x = "",
      y = "No. of Exposures",
      fill = "",
      title = "Normality of Exposure Variables",
      subtitle = "Shapiro-Wilk Test"
    )
  
  return(transform_plot)
  
}


