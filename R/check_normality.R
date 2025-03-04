#' Check Normality of Exposure Variables
#'
#' This function performs a Shapiro-Wilk test for normality on numeric exposure variables in the `colData` of an `expomicset` object. It provides a summary table and a visualization of the normality test results. Optionally, the results can be added to the metadata of `expomicset`.
#'
#' @param expomicset a `MultiAssayExperiment` object containing exposure data in `colData`.
#' @param action a character string. Can be "add" or "get". "add" will add to the `metadata` of `expomicset` and "get" will return the result
#'
#' @return If `action = "add"`, returns the modified `expomicset` with normality results in its metadata.
#' If `action = "get"`, returns a list with:
#' \item{norm_df}{A data frame containing the Shapiro-Wilk test results.}
#' \item{norm_plot}{A ggplot object visualizing the normality test results.}
#'
#' @details
#' Numeric columns in `colData` are tested using the Shapiro-Wilk test. We remove any variables that are constant, meaning they have the same values across samples before testing. Addtionally a bar plot is generated, which summarizes the proportion of normal and non-normal exposures.
#'
#' @import MultiAssayExperiment tidyverse broom ggpubr ggsci
#' @export
#'
#' @examples
#' \dontrun{
#'   # Checking the normality of the exposure variables
#'   expom <- expom |> 
#'       check_normality(expomicset, action = "get")
#' }
check_normality <- function(expomicset,
                            action="add") {
  require(tidyverse)
  require(broom)
  require(ggpubr)
  require(ggsci)
  require(MultiAssayExperiment)
  
  message("Checking Normality Using Shapiro-Wilk Test")
  
  # Extract numeric exposure data from colData
  exposure_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    dplyr::select_if(is.numeric) |>
    dplyr::select_if(function(x) !all(x == x[1]))  # Remove constant columns
  
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
    dplyr::mutate(exposure = colnames(exposure_data))
  
  # Create normality plot
  norm_plot <- table(norm_df$p.value > 0.05) |>
    as.data.frame() |>
    dplyr::mutate(Var1 = dplyr::case_when(
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
    ggpubr::theme_pubr(legend = "right") +
    ggsci::scale_fill_lancet() +
    ggsci::scale_color_lancet(guide = FALSE) +
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
  
  if(action=="add"){
    # Add normality results to expomicset metadata
    MultiAssayExperiment::metadata(expomicset)$normality <- list(
      norm_df = norm_df,
      norm_plot = norm_plot
    )
    return(expomicset)
  } else if (action=="get"){
    # Return normality results
    return(list(norm_df = norm_df, norm_plot = norm_plot))
  } else {
    stop("Invalid action. Use 'add' to add results to expomicset or 'get' to return results.")
  }
}
