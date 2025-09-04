test_that(".load_ontologies returns a list with hpo/ecto/chebi", {
  ont <- .load_ontologies(use_demo = TRUE)
  expect_type(ont, "list")
  expect_true(all(c("hpo", "ecto", "chebi") %in% names(ont)))
})


test_that(".run_categorize_ontology assigns a root_id", {
  toy <- data.frame(
    id   = c("A","B"),
    name = c("root","child"),
    ancestors = I(list(character(0), "A")),
    parents   = I(list(character(0), "A")),
    children  = I(list("B", character(0))),
    depth     = c(0,1)
  )
  df <- data.frame(selected_ontology_id = "B")
  out <- .run_categorize_ontology(df, id_col="selected_ontology_id", ontologyDF=toy)
  expect_equal(unname(out$root_id), "A")
})

test_that(".can_use_network respects option", {
  old <- getOption("tidyexposomics.allow_network")
  options(tidyexposomics.allow_network = FALSE)
  expect_false(.can_use_network())
  options(tidyexposomics.allow_network = old)
})


test_that("build_ont_annot_app returns shiny.appobj", {
  app <- build_ont_annot_app(use_demo = TRUE)
  expect_s3_class(app, "shiny.appobj")
})

