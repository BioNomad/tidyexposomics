#' Correlate Exposure Variables
#'
#' Computes correlations between numeric, categorical, and mixed exposure variables
#' using Pearson correlation, R-squared from linear models, and Cramer's V.
#' Supports filtering based on correlation strength.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure variable data.
#' @param exposure_cols An optional character vector specifying exposure variables to include. If `NULL`, all numeric exposures are used. Default is `NULL`.
#' @param threshold A numeric value specifying the minimum absolute correlation required to retain results. Default is `0.3`.
#' @param action A character string indicating whether to return results (`"get"`) or add them to metadata (`"add"`). Default is `"add"`.
#'
#' @details
#' The function extracts numeric and categorical exposure variables from `colData(expomicset)` and computes:
#' - **Numeric-numeric correlations** using Pearson correlation.
#' - **Numeric-categorical associations** using R-squared from linear models.
#' - **Categorical-categorical associations** using Cramer's V.
#'
#' Results are filtered based on the absolute correlation threshold and stored in metadata if `action = "add"`.
#'
#' @return If `action = "add"`, returns the modified `expomicset` with correlation results stored in metadata.
#' If `action = "get"`, returns a list containing:
#' \item{correlation_table}{A data frame with all computed correlations.}
#' \item{filtered_table}{A data frame containing only correlations above the threshold.}
#'
#' @examples
#' \dontrun{
#' results <- correlate_exposures(
#'   expomicset = expom,
#'   exposure_cols = c("PM2.5", "NO2", "lead"),
#'   threshold = 0.3,
#'   action = "get"
#' )
#' }
#'
#' @export
correlate_exposures <- function(
    expomicset,
    exposure_cols=NULL,
    threshold = 0.3,
    action = "add") {
  message("Characterizing exposure variables")

  # Extract and preprocess colData
  exposure_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    select_if(~ !all(. == .[1]))

  if(is.null(exposure_cols)) {
    exposure_cols <- colnames(exposure_data)
  } else {
    exposure_cols <- intersect(exposure_cols, colnames(exposure_data))
    exposure_data <- exposure_data |>
      dplyr::select(dplyr::all_of(exposure_cols))
  }

  # Identify categorical and numeric columns
  categorical_vars <- exposure_data |>
    dplyr::select_if(is.character) |>
    colnames()
  numeric_vars <- exposure_data |>
    dplyr::select_if(is.numeric) |>
    colnames()

  # Generate combinations of numeric-categorical and categorical-categorical
  num_cat_combinations <- expand.grid(
    numeric_vars,
    categorical_vars) |>
    setNames(c("var1", "var2")) |>
    dplyr::mutate(dplyr::across(
      dplyr::everything(),
      as.character)) |>
    dplyr::mutate(var1 = pmin(var1, var2),
                  var2 = pmax(var1, var2)) |>
    dplyr::filter(var1 != var2) |>
    dplyr::distinct()

  cat_cat_combinations <- expand.grid(categorical_vars, categorical_vars) |>
    setNames(c("var1", "var2")) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        as.character)) |>
    dplyr::mutate(var1 = pmin(var1, var2),
                  var2 = pmax(var1, var2)) |>
    dplyr::filter(var1 != var2) |>
    dplyr::distinct()

  message("Calculating Numeric-numeric correlations")

  # Numeric-to-numeric correlations
  num_num_corr <- cor(exposure_data |>
                        dplyr::select(
                          dplyr::all_of(numeric_vars)),
                      method = "pearson") |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(-var1,
                        names_to = "var2",
                        values_to = "correlation") |>
    dplyr::mutate(var1 = pmin(var1, var2),
                  var2 = pmax(var1, var2)) |>
    dplyr::filter(var1 != var2) |>
    dplyr::distinct()

  message("Calculating Numeric-categorical and Categorical-categorical correlations")

  # Numeric-to-categorical correlations (R-squared from lm)
  num_cat_corr <- num_cat_combinations |>
    dplyr::mutate(correlation = purrr::map2_dbl(
      var1, var2, ~ {
      num_var <- if (is.numeric(exposure_data[[.x]])) exposure_data[[.x]] else exposure_data[[.y]]
      cat_var <- if (is.character(exposure_data[[.x]])) exposure_data[[.x]] else exposure_data[[.y]]
      lm_result <- lm(num_var ~ cat_var)
      summary(lm_result)$r.squared |> replace_na(NA)
    })) |>
    dplyr::filter(!is.na(correlation))


  message("Calculating Categorical-categorical correlations")

  # Categorical-to-categorical correlations (Cramer's V)
  cat_cat_corr <- cat_cat_combinations |>
    dplyr::mutate(correlation = purrr::map2_dbl(
      var1,
      var2, ~ {
      table_data <- table(
        exposure_data[[.x]],
        exposure_data[[.y]])
      if (any(rowSums(table_data) == 0) || any(colSums(table_data) == 0)) return(NA)

      # Perform chi-squared test with simulation
      chi2 <- chisq.test(table_data, simulate.p.value = TRUE)

      # Calculate Cramer's V
      sqrt(chi2$statistic / sum(table_data) / min(nrow(table_data) - 1,
                                                  ncol(table_data) - 1))
    })) |>
    dplyr::filter(!is.na(correlation))

  # Combine all correlations
  all_corr <- dplyr::bind_rows(num_num_corr,
                               num_cat_corr,
                               cat_cat_corr) |>
    dplyr::mutate(abs_correlation = abs(correlation)) |>
    dplyr::inner_join(MultiAssayExperiment::metadata(expomicset)$var_info |>
                 as.data.frame() |>
                 dplyr::select(variable, category) |>
                 dplyr::rename(category_1 = category),
               by = c("var1" = "variable")) |>
    dplyr::inner_join(MultiAssayExperiment::metadata(expomicset)$var_info |>
                 as.data.frame() |>
                 dplyr::select(variable, category) |>
                 dplyr::rename(category_2 = category),
               by = c("var2" = "variable"))

  # Filter correlations by threshold
  filtered_corr <- all_corr |>
    dplyr::filter(abs_correlation > threshold)

  if(action=="add"){
    # Save results and heatmap in metadata
    MultiAssayExperiment::metadata(expomicset)$exposure_correlation <- list(
      correlation_table = all_corr,
      filtered_table = filtered_corr
    )

    return(expomicset)
  }else if (action=="get"){
    return(list(
      correlation_table = all_corr,
      filtered_table = filtered_corr
    ))
  }else{
    stop("Invalid action. Use 'add' or 'get'.")
  }
}
