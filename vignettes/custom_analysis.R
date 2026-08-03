## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-libraries,message=FALSE,warning=FALSE-------------------------------
# Load Libraries
library(tidyverse)
library(tidyexposomics)

## ----load-data----------------------------------------------------------------
# Load example data
data("tidyexposomics_example")

# Create exposomic set object
expom <- create_exposomicset(
    codebook = tidyexposomics_example$annotated_cb,
    exposure = tidyexposomics_example$meta,
    omics = list(
        "Gene Expression" = tidyexposomics_example$exp_filt,
        "Methylation" = tidyexposomics_example$methyl_filt
    ),
    row_data = list(
        "Gene Expression" = tidyexposomics_example$exp_fdata,
        "Methylation" = tidyexposomics_example$methyl_fdata
    )
)

## ----pivot-sample-------------------------------------------------------------

# Pivot sample data to a tibble
expom |>
    pivot_sample() |>
    head()

## ----pivot-sample-example-----------------------------------------------------

# Count the number of asthma cases by sex status
expom |>
    pivot_sample() |>
    group_by(hs_asthma, e3_sex_None) |>
    reframe(n = n())

## ----pivot-feature------------------------------------------------------------

# Pivot feature data to a tibble
expom |>
    pivot_feature() |>
    head()

## ----pivot-feature-example----------------------------------------------------

# Count the number of features per omic layer
expom |>
    pivot_feature() |>
    group_by(.exp_name) |>
    summarise(n = n())

## ----pivot-exp----------------------------------------------------------------

# Pivot experiment data to a tibble
expom |>
    pivot_exp(
        exp_name = "Gene Expression",
        features = "TC01004453.hg.1"
    ) |>
    head()

## ----pivot-exp-example, fig.height=4, fig.width=4, fig.align='center',fig.cap="IL23R gene expression levels (probe TC01004453.hg.1) by asthma status."----

# Plot expression of IL23R
expom |>
    pivot_exp(
        exp_name = "Gene Expression",
        features = "TC01004453.hg.1"
    ) |>
    ggplot(aes(
        x = as.character(hs_asthma),
        y = log2(counts + 1),
        color = as.character(hs_asthma),
        fill = as.character(hs_asthma)
    )) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(alpha = 0.1) +
    ggpubr::geom_pwc(label = "{p.adj.format}{p.adj.signif}") +
    theme_minimal() +
    ggpubr::rotate_x_text(angle = 45) +
    ggsci::scale_color_cosmic() +
    ggsci::scale_fill_cosmic() +
    labs(
        x = "",
        y = expression(Log[2] * "Abd."),
        fill = "Asthma Status",
        color = "Asthma Status"
    )

## ----session_info,collapse = TRUE---------------------------------------------

sessionInfo()

