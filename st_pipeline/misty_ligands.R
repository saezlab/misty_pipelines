# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing of breast cancer visium slides
#' It normalizes the data using SCT transform,
#' calculates TF and PROGENy activities and
#' defines a matrix of expressed ligands
#' 

library(cowplot)
source("slide_processing.R") #Check dependencies

breast_A_1 = process_visium(dir_path = "breast_A_1")
saveRDS(breast_A_1, file = "breast_A_1.rds")
breast_A_2 = process_visium(dir_path = "breast_A_2")
saveRDS(breast_A_2, file = "breast_A_2.rds")


# General features from slides
# Mean genes
mean(c(mean(breast_A_1@meta.data$nFeature_Spatial),
       mean(breast_A_2@meta.data$nFeature_Spatial)))

mean(c(mean(breast_A_1@meta.data$nFeature_progeny),
mean(breast_A_2@meta.data$nFeature_progeny)))

mean(c(ncol(breast_A_1),ncol(breast_A_2)))

mean(c(mean(breast_A_1@meta.data$nFeature_ligands),
       mean(breast_A_2@meta.data$nFeature_ligands)))

# Intersection of expressed ligands
ligands = intersect(rownames(breast_A_1@assays$ligands@data),
                    rownames(breast_A_2@assays$ligands@data))

# Run multiple l parameters

ls = c(2,5,10,20,50,100)

test_A1 = lapply(ls, run_ligand_MISTy, seurat_visium_obj = breast_A_1,
                 ligs_features = ligands, out_dir_name = "mrun_breast1")

get_optimal(out_dir_name = "mrun_breast1",ls = ls)

test_A2 = lapply(ls, run_ligand_MISTy, seurat_visium_obj = breast_A_2,
                 ligs_features = ligands, out_dir_name = "mrun_breast2")

get_optimal(out_dir_name = "mrun_breast2",ls = ls)

# Generate output for paper
results_folder = c("mrun_breast1_optim",
                   "mrun_breast2_optim")

MISTy_out = MISTy_aggregator(results_folder = results_folder,
                 p.cutoff = 0.05)

# In average which view contributed the most?

MISTy_out$coefs %>% group_by(view) %>%
  summarize(mean(value))

MISTy_out$coefs %>% group_by(view) %>%
  filter(value == max(value))

MISTy_out$coefs %>% group_by(view) %>%
  arrange(desc(value)) %>% slice(1:3)
