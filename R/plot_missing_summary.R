plot_missing_summary <- function(
    expomicset,
    threshold=20){
  require(tidyverse)
  require(ggpubr)
  require(ggsci)
  
  # Plot missing data summary
  qc_res <- expomicset@metadata$na_qc
  
  # loop through each omic and determine which layers have missing data
  missing_data <- imap(expomicset@metadata$na_qc, ~ {
    #df <- .x$vars_to_exclude_omics_sum
    df <- .x$all_var_sum |> 
      filter(pct_miss>threshold)
    
    # Ensure df is a data frame
    if (is.null(df)) {
      return(NULL)  # Skip this entry if NULL
    }
    
    if (nrow(df) > 0) {
      df <- df |> mutate(exp_name = .y)
    }
    
    return(df)  # Always return df, even if empty
  }) |>
    bind_rows() |> 
    group_by(exp_name) |>
    summarise(missingness = n()) 
  
  # add in other assay names and say their missingness is 0
  exp_names <- expomicset@metadata$na_qc |> names()
  missing_data <- missing_data |> 
    bind_rows(
      data.frame(
        exp_name = exp_names[!exp_names %in% missing_data$exp_name],
        missingness = 0)) |> 
    mutate(exp_name=case_when(
      exp_name == "exposure" ~ "Exposure",
      .default = exp_name
    )) |> 
    mutate(exp_name = paste(exp_name, " (", missingness, ")", sep = ""))
  
  # plot missing data summary
  missing_data |> 
    ggplot(aes(y = fct_reorder(exp_name, missingness), 
               x = missingness, 
               fill = exp_name)) +
    geom_bar(stat = "identity",alpha=0.7) +
    geom_segment(aes(
      x = missingness,                    
      xend = missingness,                    
      y = as.numeric(fct_reorder(exp_name, missingness)) - 0.45,
      yend = as.numeric(fct_reorder(exp_name, missingness)) + 0.45,
      color = exp_name,
    ), size = 1) +
    scale_fill_npg()+
    scale_color_npg()+
    labs(title = paste("No. of Features Over ", threshold, "%"," Threshold",sep=""),
         y = "",
         x = "No. of Features") +
    theme_pubr(legend="none") 
  
  #return(missing_data)
}