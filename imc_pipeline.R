library(future)
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(MISTy)

data <- list.dirs("/data/IMC_breastcancer", recursive = FALSE)
plan(multiprocess, workers = 4)

threshold <- 11.236


ls <- c(25, 50, 100, 200, 400)

ls %>% walk(function(l) {
  data %>% walk(function(d) {
    expr <- read_delim(paste0(d, "/expressions.txt"), delim = " ", col_types = cols())
    colnames(expr) <- make.names(colnames(expr))
    pos <- read_csv(paste0(d, "/positions.txt"), col_names = c("x", "y"), col_types = cols())

    views <- create_initial_view(expr) %>%
      add_juxtaview(pos, threshold) %>%
      add_paraview(pos, l^2)

    run_misty(views, results.folder = paste0(
      "results/imc_bc_", l^2, "/",
      str_extract(d, "[ABC][a-zA-Z0-9]+"), "/")
    )
  })
})
