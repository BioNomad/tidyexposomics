# 
# pacman::p_load(rvest,xml2,tidyverse)
# 
# search <- "arsenic"
# # URL of the page
# url <- paste0(
#   "https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search=",
#   search,
#   "&display=200000&zsort=samples")
# 
# # Read the webpage
# webpage <- read_html(url)
# 
# # Extract the table
# table_node <- html_node(webpage, "table")  # Find the table element
# table_df <- html_table(table_node, fill = TRUE)  # Convert to a dataframe
# 
# table_df <- table_df |> 
#   mutate(`Series type(s)`=gsub(
#     "\n                    \n                    \n                      ",
#     ",",
#     `Series type(s)`)) |> 
#   mutate(Supplementary=gsub(
#     " \n                    \n                    \n                      \n                        \n                      ",
#     ",",
#     Supplementary)) 
# 


scrape_geo_data <- function(search_term) {
  require(rvest)
  require(xml2)
  require(tidyverse)
  
  # Construct search URL
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search=",
    search_term,
    "&display=200000&zsort=samples"
  )
  
  # Read the webpage
  webpage <- read_html(url)
  
  # Extract the table
  table_node <- html_node(webpage, "table")  
  table_df <- html_table(table_node, fill = TRUE)  
  
  # Clean 'Series type(s)' column
  table_df <- table_df %>% 
    mutate(`Series type(s)` = str_replace_all(
      `Series type(s)`, "\\s{2,}", ","))  # Replace excessive spaces with commas
  
  # Clean 'Supplementary' column and split into new columns
  table_df <- table_df %>% 
    mutate(Supplementary = str_replace_all(
      Supplementary, "\\s{2,}", ","))  # Replace excessive spaces with commas
  
  return(table_df)
}

# a=scrape_geo_data("autism") |> 
#   filter(`Organism(s)` == "Homo sapiens")
