library(future)
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(mistyR)

data <- list.dirs("../data/imc_small_breastcancer/", recursive = FALSE)
plan(multiprocess, workers = 4)

threshold <- 11.236

perm <- seq(0,10)

ls <- c(25, 50, 100, 200, 400)

ls %>% walk(function(l) {
  data %>% walk(function(d) {
    expr <- read_delim(paste0(d, "/expressions.txt"), delim = " ", col_types = cols())
    colnames(expr) <- make.names(colnames(expr))
    pos <- read_csv(paste0(d, "/positions.txt"), col_names = c("x", "y"), col_types = cols())
    
    perm %>% walk(function(p){
      
      if(p != 0){
        set.seed(p)
        perm.pos <- pos[sample(nrow(pos)), ]
      } else {
        perm.pos <- pos
      }
        
      views <- create_initial_view(expr) %>%
        add_juxtaview(perm.pos, threshold, cached = FALSE) %>%
        add_paraview(perm.pos, l, cached = FALSE)
  
      run_misty(views, results.folder = paste0(
        "../results/imc_small_perm",p,"/imc_small_", l, "/",
        str_extract(d, "[ABC][a-zA-Z0-9]+"), "/"),
        cached = FALSE
      )
    })
  })
})
