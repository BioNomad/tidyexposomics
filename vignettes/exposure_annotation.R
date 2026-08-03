## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----ontology-app,eval=FALSE--------------------------------------------------
# # Launch the shiny app to annotate exposure variables
# shiny::runApp(build_ont_annot_app())

## ----ontology-app-img, echo=FALSE,out.width="150%",fig.align='center',fig.cap="Screenshot of the ontology annotation app. The sidebar has an upload button to load your exposure metadata file. The main panel displays the uploaded exposure metadata where users can select variables to annotate and categorize. The sidebar also contains buttons to apply annotations/categorizations and download the annotated file."----
knitr::include_graphics("./ont_annot.png")

## ----session_info,collapse = TRUE---------------------------------------------
sessionInfo()

