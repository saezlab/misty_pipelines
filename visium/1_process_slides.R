# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run the basic Seurat pipeline for spatial transcriptomics
#' folder
#' |
#' sample---
#'          |
#'          ---spatial
#'          ---filtered_feature_bc_matrix.h5
#' 1) Read 10x files
#' 2) Normalization with cpm and SCT
#' 3) PROGENy

library(tidyverse)
library(Seurat)
library(progeny)

#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param top: number of genes used in footprint
#' @return A Seurat object with pathway activities in progeny assay
add_path_activities <- function(visium_slide, 
                                species = "human",
                                top = 500, 
                                assay = "RNA"){
  
  if(species == "mouse"){
    model <- progeny::getModel(organism = "Mouse", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- (t(model)[, common_genes]) %*% (GetAssayData(visium_slide, assay = assay)[common_genes, ])
    # Here we scale by feature
    progeny_scores <- scale(t(as.matrix(progeny_scores)))
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
  
  }else if(species == "human"){
    
    model <- progeny::getModel(organism = "Human", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- (t(model)[, common_genes]) %*% (GetAssayData(visium_slide, assay = assay)[common_genes, ])
    # Here we scale by feature
    progeny_scores <- scale(t(as.matrix(progeny_scores)))
    
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
  }
  
  return(visium_slide)
}


## Process all data

param_df <- tibble("sample_name" = list.dirs("./data",full.names = F, recursive = F),
                   "data_dir" = list.dirs("./data",full.names = T, recursive = F)) %>%
  dplyr::mutate("out_folder" = paste0("./processed_visium/", sample_name, ".rds"))

walk2(param_df$data_dir, param_df$out_folder, function(dir_name, out_dir_name) {
  
  print(dir_name)
  
  sample_name <- gsub("./data/", "", dir_name)
  
  visium_slide = Load10X_Spatial(data.dir = dir_name) %>%
    SCTransform(. , assay = "Spatial", verbose = F) %>%
    NormalizeData(., 
                  normalization.method = 'LogNormalize', 
                  scale.factor = 10000, 
                  verbose = FALSE) %>%
    add_path_activities(.,
                        species = "human",
                        top = 1000, 
                        assay = "SCT")
  
  saveRDS(visium_slide, file = out_dir_name)
  
  DefaultAssay(visium_slide) <- "progeny"
  
  
  pdf(paste0("./results/progeny_spatial_plts/", sample_name, ".pdf"), height = 12, width = 10)
  
  all_paths <- SpatialFeaturePlot(visium_slide, 
                                  features = rownames(visium_slide), 
                                  ncol = 4,
                                  stroke = 0,
                                  max.cutoff = "q99")
  plot(all_paths)
  
  dev.off()
  
})

