#' Assess Normality of Exposure Variables
#'
#' Performs Shapiro-Wilk tests to check the normality of numeric exposure variables in `colData`
#' of a `MultiAssayExperiment` object.
#'
#' @param expomicset A `MultiAssayExperiment` object containing exposure data in `colData`.
#' @param action A character string specifying whether to store (`"add"`) or return (`"get"`) the results. Default is `"add"`.
#'
#' @details
#' This function:
#' - Extracts **numeric, non-constant** exposure variables from `colData`.
#' - Runs **Shapiro-Wilk tests** to assess normality.
#' - Summarizes the number of normally and non-normally distributed exposures.
#' - Generates a bar plot visualizing the normality results.
#' - **Output Handling**:
#'   - `"add"`: Stores results in `metadata(expomicset)$normality`.
#'   - `"get"`: Returns a list containing the normality test results and plot.
#'
#' @return A `MultiAssayExperiment` object with normality results added to metadata (if `action = "add"`) or a list with:
#' \item{norm_df}{A data frame of Shapiro-Wilk test results for each exposure variable.}
#' \item{norm_plot}{A ggplot object showing the distribution of normal and non-normal exposures.}
#'
#' @examples
#' \dontrun{
#' expom <- run_normality_check(expomicset = expom, action = "add")
#' norm_results <- run_normality_check(expomicset = expom, action = "get")
#' }
#'
#' @export
run_normality_check <- function(expomicset,
                            action="add") {
  require(ggplot2)

  message("Checking Normality Using Shapiro-Wilk Test")

  # Extract numeric exposure data from colData
  exposure_data <- MultiAssayExperiment::colData(expomicset) |>
    as.data.frame() |>
    dplyr::select_if(is.numeric) |>
    dplyr::select_if(function(x) !all(x == x[1]))

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

  # Summarize normality results
  norm_summary <- norm_df |>
    dplyr::summarise(
      "Normal" = sum(p.value > 0.05),
      "Not Normal" = sum(p.value <= 0.05)
    ) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("var") |>
    setNames(c("var","value"))

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
      norm_summary = norm_summary
    )
    return(expomicset)
  } else if (action=="get"){
    # Return normality results
    return(list(norm_df = norm_df, norm_plot = norm_plot))
  } else {
    stop("Invalid action. Use 'add' to add results to expomicset or 'get' to return results.")
  }
}
