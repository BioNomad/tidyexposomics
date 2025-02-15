check_normality <- function(expOmicSet) {
  require(tidyverse)
  require(broom)
  require(ggpubr)
  require(ggsci)
  
  message("Checking Normality Using Shapiro-Wilk Test")
  
  # Extract numeric exposure data from colData
  exposure_data <- colData(expOmicSet) |>
    as.data.frame() |>
    select_if(is.numeric) |>
    select_if(function(x) !all(x == x[1]))  # Remove constant columns
  
  if (ncol(exposure_data) == 0) {
    stop("No numeric or non-constant exposure variables found for normality testing.")
  }
  
  # Perform Shapiro-Wilk test for normality
  norm_df <- exposure_data |>
    apply(2, function(x) {
      shapiro.test(x) |>
        broom::tidy()
    }) |>
    (\(x) do.call(rbind, x))() |>
    mutate(exposure = colnames(exposure_data))
  
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
  
  # Print normality summary
  num_normal <- table(norm_df$p.value > 0.05)["TRUE"] |> as.numeric()
  num_not_normal <- table(norm_df$p.value < 0.05)["TRUE"] |> as.numeric()
  
  message(ifelse(is.na(num_normal), 0, num_normal),
          " Exposure Variables are Normally Distributed")
  
  message(ifelse(is.na(num_not_normal), 0, num_not_normal),
          " Exposure Variables are NOT Normally Distributed")
  
  # Save results in expOmicSet
  metadata(expOmicSet)$normality <- list(
    norm_df = norm_df,
    norm_plot = norm_plot
  )
  
  return(expOmicSet)
}
