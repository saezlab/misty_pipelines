# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Building a MISTy pipeline to predict pathway
#' activities using ligands from Omnipath
#' 
#' Generates results from permuted data to compare
#' performances
#' 

library(OmnipathR)
library(tidyverse)
library(Seurat)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

breast_A_1 = readRDS("breast_A_1.rds")
breast_A_2 = readRDS("breast_A_2.rds")

# Importing ligand-receptor network
lig_rec = import_intercell_network(
  interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
  transmitter_param = list(parent = 'ligand'),
  receiver_param = list(parent = 'receptor')
)

# Selecting important features
lig_rec_red = lig_rec %>% dplyr::ungroup() %>%
  dplyr::select(target_genesymbol,source_genesymbol,
                is_directed,category_intercell_target)

# Making a set of unique ligands and receptors
ligands = unique(lig_rec_red$source_genesymbol)

#Checking which ligand receptors are in both slides (Coverage 30%)
all_slide_markers = map(list(breast_A_1,breast_A_2), function(visium_slide){
  slide_markers = ligands[ligands %in% rownames(visium_slide)]
  gex = as.matrix(visium_slide@assays$SCT@data)[slide_markers,]
  coverage = rowSums(gex>0)/ncol(gex)
  slide_markers = names(coverage[coverage>=.05])
  return(slide_markers)
})

all_slide_markers = intersect(all_slide_markers[[1]],
                              all_slide_markers[[2]])

#Defining constant parameter of the pipeline
view_assays = list("main" = "progeny",
                   "para_path" = "progeny",
                   "para_lig" = "SCT")

view_features = list("main" = NULL, #uses all
                     "para_path" = NULL, #uses all
                     "para_lig" = unique(all_slide_markers))

view_types = list("main" = "intra", #uses all
                  "para_path" = "para", #uses all
                  "para_lig" = "para")

view_params = list("main" = NULL, #uses all
                   "para_path" = NULL, #uses all
                   "para_lig" = NULL)

ls = c(2,5,10)

# Generate pipeline

# plan(multiprocess, workers = workers)
# clear_cache()

run_misty_pipeline = function(visium_slide, ls, out_prealias = "./revisions_misty/results/A1model"){
  
  for(l in ls){
    
    misty_out = sprintf(paste0(out_prealias,"_%s"), l^2)
    
    view_params = list("main" = NULL, #uses all
                       "para_path" = l, #uses all
                       "para_lig" = l)
    
    MISTy_seurat(visium_slide = visium_slide,
                 view_assays = view_assays,
                 view_features = view_features,
                 view_types = view_types,
                 view_params = view_params,
                 spot_ids = NULL,
                 out_alias = misty_out,
                 workers = 4)
  }
  
  get_optimal(out_dir_name = out_prealias,ls = ls)

}

# Run original pipelines
# A1

run_misty_pipeline(visium_slide = breast_A_1, 
                   ls = ls, 
                   out_prealias = "./revisions_misty/results/A1model")

for(i in 1:5){
  print(i)
  set.seed(i)
  
  visium_rndm = breast_A_1
  rndm_geometry = breast_A_1@images$slice1@coordinates
  # Here we shuffle the geometry
  rownames(rndm_geometry) = sample(rownames(rndm_geometry))
  visium_rndm@images$slice1@coordinates = rndm_geometry
  
  run_misty_pipeline(visium_slide = visium_rndm, 
                     ls = ls, 
                     out_prealias = sprintf("./revisions_misty/results/A1model_rndm_%s",i))
  
  
}


# A2
run_misty_pipeline(visium_slide = breast_A_2, 
                   ls = ls, 
                   out_prealias = "./revisions_misty/results/A2model")

for(i in 1:5){
  print(i)
  set.seed(i)
  
  visium_rndm = breast_A_2
  rndm_geometry = breast_A_2@images$slice1@coordinates
  # Here we shuffle the geometry
  rownames(rndm_geometry) = sample(rownames(rndm_geometry))
  visium_rndm@images$slice1@coordinates = rndm_geometry
  
  run_misty_pipeline(visium_slide = visium_rndm, 
                     ls = ls, 
                     out_prealias = sprintf("./revisions_misty/results/A2model_rndm_%s",i))
  
  
}