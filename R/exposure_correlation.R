exposure_correlation <- function(
    expomicset, 
    exposure_cols=NULL,
    threshold = 0.3,
    action = "add") {
  library(tidyverse)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  
  message("Characterizing exposure variables")
  
  # Extract and preprocess colData
  exposure_data <- colData(expomicset) |>
    as.data.frame() |>
    select_if(~ !all(. == .[1]))  # Keep columns with more than one unique value
  
  if(is.null(exposure_cols)) {
    exposure_cols <- colnames(exposure_data)
  } else {
    exposure_cols <- intersect(exposure_cols, colnames(exposure_data))
    exposure_data <- exposure_data |>
      select(all_of(exposure_cols))
  }
  
  # Identify categorical and numeric columns
  categorical_vars <- exposure_data |>
    select_if(is.character) |>
    colnames()
  numeric_vars <- exposure_data |>
    select_if(is.numeric) |>
    colnames()
  
  # Generate combinations of numeric-categorical and categorical-categorical
  num_cat_combinations <- expand.grid(numeric_vars, categorical_vars) |>
    setNames(c("var1", "var2")) |>
    mutate(across(everything(), as.character)) |>
    mutate(var1 = pmin(var1, var2), var2 = pmax(var1, var2)) |>
    filter(var1 != var2) |>
    distinct()
  
  cat_cat_combinations <- expand.grid(categorical_vars, categorical_vars) |>
    setNames(c("var1", "var2")) |>
    mutate(across(everything(), as.character)) |>
    mutate(var1 = pmin(var1, var2), var2 = pmax(var1, var2)) |>
    filter(var1 != var2) |>
    distinct()
  
  message("Calculating Numeric-numeric correlations")
  
  # Numeric-to-numeric correlations
  num_num_corr <- cor(exposure_data |>
                        dplyr::select(all_of(numeric_vars)),
                      method = "pearson") |>
    as.data.frame() |>
    rownames_to_column("var1") |>
    pivot_longer(-var1, names_to = "var2", values_to = "correlation") |>
    mutate(var1 = pmin(var1, var2), var2 = pmax(var1, var2)) |>
    filter(var1 != var2) |>
    distinct()
  
  message("Calculating Numeric-categorical and Categorical-categorical correlations")
  # Numeric-to-categorical correlations (R-squared from lm)
  num_cat_corr <- num_cat_combinations |>
    mutate(correlation = map2_dbl(var1, var2, ~ {
      num_var <- if (is.numeric(exposure_data[[.x]])) exposure_data[[.x]] else exposure_data[[.y]]
      cat_var <- if (is.character(exposure_data[[.x]])) exposure_data[[.x]] else exposure_data[[.y]]
      lm_result <- lm(num_var ~ cat_var)
      summary(lm_result)$r.squared |> replace_na(NA)
    })) |>
    filter(!is.na(correlation))
  
  
  message("Calculating Categorical-categorical correlations")
  # Categorical-to-categorical correlations (Cramer's V)
  cat_cat_corr <- cat_cat_combinations |>
    mutate(correlation = map2_dbl(var1, var2, ~ {
      table_data <- table(exposure_data[[.x]], exposure_data[[.y]])
      if (any(rowSums(table_data) == 0) || any(colSums(table_data) == 0)) return(NA)
      chi2 <- chisq.test(table_data, simulate.p.value = TRUE)
      sqrt(chi2$statistic / sum(table_data) / min(nrow(table_data) - 1, ncol(table_data) - 1))
    })) |>
    filter(!is.na(correlation))
  
  # Combine all correlations
  all_corr <- bind_rows(num_num_corr, num_cat_corr, cat_cat_corr) |>
    mutate(abs_correlation = abs(correlation)) |>
    inner_join(expomicset@metadata$var_info |>
                 as.data.frame() |>
                 dplyr::select(variable, category) |>
                 dplyr::rename(category_1 = category),
               by = c("var1" = "variable")) |>
    inner_join(expomicset@metadata$var_info |>
                 as.data.frame() |>
                 dplyr::select(variable, category) |>
                 dplyr::rename(category_2 = category),
               by = c("var2" = "variable"))
  
  # Filter correlations by threshold
  filtered_corr <- all_corr |>
    filter(abs_correlation > threshold)
  
  # Generate heatmap
  heatmap <- ggplot(filtered_corr, aes(x = var1, y = var2, fill = abs_correlation)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "magenta4") +
    theme_void() +
    rotate_x_text(angle = 90) +
    theme(
      strip.text.x = element_text(angle = 90),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    facet_grid(category_2 ~ category_1, scales = "free", space = "free") +
    labs(fill = "Abs.(Correlation)")
  
  if(action=="add"){
    # Save results and heatmap in metadata
    metadata(expomicset)$exposure_correlation <- list(
      correlation_table = all_corr,
      filtered_table = filtered_corr,
      heatmap = heatmap
    )
    
    return(expomicset)
  }else if (action=="get"){
    return(list(
      correlation_table = all_corr,
      filtered_table = filtered_corr,
      heatmap = heatmap
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}



# # select data with more than one value type
# a <- expom_7 |> 
#   colData() |> 
#   as.data.frame() |>
#   select_if( ~ !all(.==.[1]))
# 
# # grab categorical columns
# b <- a |> 
#   select_if(is.character) |>
#   colnames()
# 
# # grab numeric columns
# c <- a |> 
#   select_if(is.numeric) |>
#   colnames()
# 
# # combinations of numeric-categorical
# d <- expand.grid(c, b) |>
#   setNames(c("var1", "var2")) |>
#   mutate(across(everything(), as.character)) |> 
#   mutate(
#     var1 = pmin(var1, var2),  
#     var2 = pmax(var1, var2)   
#   ) |> 
#   filter(var1 != var2) |> 
#   distinct()  
# 
# # combinations of categorical-categorical
# e <- expand.grid(b,b) |>
#   setNames(c("var1", "var2")) |> 
#   mutate(across(everything(), as.character)) |> 
#   mutate(
#     var1 = pmin(var1, var2),  
#     var2 = pmax(var1, var2)   
#   ) |> 
#   filter(var1 != var2) |> 
#   distinct()  
# 
# # correlation between numeric data using Hsmisc
# f <- cor(a |> 
#            select(c), 
#          method = "pearson") |>
#   as.data.frame() |> 
#   rownames_to_column("var1") |> 
#   pivot_longer(-var1,
#                names_to = "var2",
#                values_to = "correlation") |> 
#   mutate(
#     var1 = pmin(var1, var2),  
#     var2 = pmax(var1, var2)   
#   ) |> 
#   filter(var1 != var2) |> 
#   distinct()  
# 
# # correlation between categorical and numeric using an lm r^2
# g <- d |> 
#   mutate(correlation = map2_dbl(var1, var2, ~{
#     
#     # Always fit the numeric variable as the dependent variable
#     num_var <- if (is.numeric(a[[.x]])) a[[.x]] else a[[.y]]
#     cat_var <- if (is.factor(a[[.x]])) a[[.x]] else a[[.y]]
#     
#     # grab the R^2 value
#     lm_result <- lm(num_var ~ cat_var)
#     r_squared <- summary(lm_result)$r.squared
#     
#     if (is.na(r_squared) || r_squared == 1) {
#       return(NA)
#     } else {
#       return(r_squared)
#     }
#   })) |> 
#   filter(!is.na(correlation))
# 
# # correlation between categorical and categorical using Cramer's V
# h <- e |> 
#   mutate(correlation = map2_dbl(var1, var2, ~{
#     table_data <- table(a[[.x]], a[[.y]])
#     chi2 <- chisq.test(table_data, simulate.p.value = TRUE)
#       return(sqrt(chi2$statistic / sum(table_data) / min(nrow(table_data) - 1, ncol(table_data) - 1)))
#   })) |> 
#   filter(!is.na(correlation))
# 
# # combine all correlations
# i <- bind_rows(f, g, h) |> 
#   #filter(abs(correlation) > 0.3) |> 
#   mutate(var1 = factor(var1),
#          var2 = factor(var2)) |> 
#   mutate(abs_correlation=abs(correlation)) |> 
#   inner_join(expom_7@metadata$var_info |> 
#                as.data.frame() |> 
#                dplyr::select(variable,category) |> 
#                dplyr::rename(category_1=category),
#              by=c("var1"="variable")) |> 
#   inner_join(expom_7@metadata$var_info |> 
#                as.data.frame() |> 
#                dplyr::select(variable,category) |> 
#                dplyr::rename(category_2=category),
#              by=c("var2"="variable")) 
# 
# i |> 
#   ggplot(aes(x = var1, 
#              y = var2, 
#              fill = abs_correlation)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "magenta4") +
#   theme_void() +
#   rotate_x_text(angle = 90) +
#   theme(
#     strip.text.x = element_text(angle = 90),
#     axis.text.x = element_blank(),axis.text.y = element_blank()
#   ) +
#   facet_grid(category_2 ~ category_1, 
#              scales = "free", 
#              space = "free") +
#   labs(
#     fill = "Abs.(Correlation)"
#   )
# 
# # use tidyheatmap to visualize
# # come up with annotation colors for unique(category_2)
# # come up with annotation colors for unique(category_1)
# # 
# # library(tidyheatmaps)
# # tidyheatmap(i,
# #             rows = var1,
# #             columns = var2,
# #             values = abs_correlation,
# #             scale = "none",
# #             colors = c("white", "red"),
# #             color_na = "white",
# #             annotation_col = c(category_2),
# #             annotation_names_row = TRUE,
# #             annotation_colors = list(category_2=get_palette("npg",k = length(unique(c(i$category_1,i$category_2)) )) |> setNames(unique(c(i$category_1,i$category_2))) 
# #             )
# # )
# 
# 
