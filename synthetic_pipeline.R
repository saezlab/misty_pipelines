library(future)
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(MISTy)

data <- list.dirs("data/synthetic/", recursive = FALSE)
plan(multiprocess, workers = 4)

l <- 10

data %>% walk(function(d) {
  all <- read_csv(paste0(d, "/random1_position_expression_real_cells.csv")) %>%
    select(-type)
  expr <- all %>% select(-row, -col)
  pos <- all %>% select(row, col)

  views <- create_initial_view(expr) %>% add_paraview(pos, l^2)

  run_misty(views, results.folder = paste0(
    "results/synthetic/",
    str_extract(d, "synthetic[0-9]+"), "/"
  ))
})
