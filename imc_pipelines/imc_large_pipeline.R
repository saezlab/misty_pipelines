library(tidyverse)
library(MISTy)
library(future)

plan(multiprocess, workers = 4)

paths <- list.files("../data/imc_large_breastcancer/single_images/", full.names = TRUE)

panel <- read_csv("../data/imc_large_breastcancer/Basel_Zuri_StainingPanel.csv") %>% 
  select(Target, FullStack) %>% filter(complete.cases(.)) %>% distinct()

neighbor.thr <- 12
l <- 100

paths %>% walk(function(path){
  pid <- str_extract(path, "[0-9]+")
  message(paste0("Processing ", pid))
  
  data <- read_csv(path)
  if(nrow(data) < 1000) return()
  
  expr <- data %>% select(contains("Intensity_MeanIntensity_FullStack"))
  
  # translate channel to marker
  markers <- tibble(channel = colnames(expr) %>% str_extract("[0-9]+") %>% as.numeric) %>% 
    left_join(panel %>% slice(-c(41, 46, 47)), by = c("channel" = "FullStack"))
  
  # cleanup both markers and expr
  to.remove <- which(markers$channel < 9 | is.na(markers$Target) | 
                       markers$channel > 47 | markers$channel %in% c(12, 25, 28, 38))
  expr <- expr %>% select(-all_of(to.remove))
  colnames(expr) <- markers %>% slice(-to.remove) %>% pull(Target) %>% make.names
  
  # get geometry
  pos <- data %>% select(Location_Center_X, Location_Center_Y)
  
  # filter nans
  pos.complete <- which(complete.cases(pos))
  
  expr <- expr %>% slice(pos.complete)
  pos <- pos %>% slice(pos.complete)
  
  #misty
  views <- create_initial_view(expr) %>%
    add_juxtaview(pos, neighbor.thr, cached = FALSE) %>%
    add_paraview(pos, l^2, cached = FALSE)
  
  run_misty(views, results.folder = paste0("../results/imc_large/",pid), cached = FALSE)
  
})

