# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de

#' Interpretation of the original model

library(OmnipathR)
library(tidyverse)
library(Seurat)
library(progeny)
library(cowplot)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

breast_A_1 <- readRDS("breast_A_1.rds")
breast_A_2 <- readRDS("breast_A_2.rds")

# Importing ligand-receptor network --------------------------------------------
lig_rec <- import_intercell_network(
  interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
  transmitter_param = list(parent = 'ligand'),
  receiver_param = list(parent = 'receptor')
  )

# Selecting important features --------------------------------------------
lig_rec_red <- lig_rec %>% dplyr::ungroup() %>%
  dplyr::select(target_genesymbol,source_genesymbol,
                is_directed,category_intercell_target)


# Optimized aggregated performance --------------------------------------------
original_run <- MISTy::collect_results(c("./revisions_misty/results/A1model_optim",
                                        "./revisions_misty/results/A2model_optim"))

# Improvement plot --------------------------------------------
measure <- "gain.R2"
misty_results <- original_run

plot_data <- misty_results$improvements.stats %>% 
  dplyr::filter(measure == !!measure)

improvement_order <- levels(reorder(plot_data$target,
                                    -plot_data$mean))

set2_orange <- "#FC8D62"

imp_results_plot <- ggplot2::ggplot(plot_data, 
                                ggplot2::aes(x = reorder(target,-mean), 
                                             y = mean)) + 
  ggplot2::geom_pointrange(ggplot2::aes(ymin = mean - sd, ymax = mean + sd)) + 
  ggplot2::geom_point(color = set2_orange) + 
  ggplot2::theme_classic() + 
  ggplot2::ylab(measure) + 
  ggplot2::xlab("Target") + 
  ggplot2::theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank())

# Contributions --------------------------------------------
plot_data <- misty_results$contributions.stats
cont_results_plot <- ggplot2::ggplot(plot_data, 
                                     ggplot2::aes(x = factor(target,
                                                             levels = improvement_order),
                                                  y = fraction)) + 
  ggplot2::geom_col(ggplot2::aes(group = view, fill = view)) + 
  ggplot2::scale_fill_brewer(palette = "Set2") + 
  ggplot2::ylab("Contribution") + 
  ggplot2::xlab("Target") + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     hjust = 1,
                                                     vjust = 0.5),
                 legend.position = "bottom")

panel_b = cowplot::plot_grid(imp_results_plot, 
                             cont_results_plot, 
                             ncol = 1, 
                             align = "v", 
                             rel_heights = c(0.5, 1))

pdf(file = "./revisions_misty/figures/panel_b.pdf", height = 6, width = 3)
plot(panel_b)
dev.off()

# Significance of change in R2 --------------------------------------------
# Adjusted p-values per image ---------------------------------------------
Rpvals_corrected <- original_run$improvements %>%
  dplyr::filter(measure == "p.R2") %>% 
  group_by(image) %>%
  dplyr::mutate(value = p.adjust(value)) %>% 
  ungroup()

RMSEpvals_corrected <- original_run$improvements %>%
  dplyr::filter(measure == "p.RMSE") %>% 
  group_by(image) %>%
  dplyr::mutate(value = p.adjust(value)) %>% 
  ungroup()

original_run$improvements <- original_run$improvements %>%
  dplyr::filter(measure != "p.RMSE",
                measure != "p.R2") %>%
  bind_rows(Rpvals_corrected,RMSEpvals_corrected)

get_misty_bic(MISTyout = original_run,
              p_value_trsh = 0.1,
              R2_trsh = 0)

# Mean contribution of each view
original_run$contributions.stats %>%
  group_by(view) %>%
  summarise(mean_contribution = mean(fraction))

# Quick summaries
improvement_summary <- original_run$improvements.stats %>%
  group_by(measure) %>%
  arrange(-mean) %>%
  nest()

view_contributions <- original_run$contributions.stats %>%
  group_by(view) %>%
  arrange(-mean) %>%
  nest()

# Plot importances of first 2 views -----------------------------------------------------------
predicted_paths <- original_run$improvements$target %>% unique()
original_run_paths <- original_run
original_run_paths$importances.aggregated <- original_run_paths$importances.aggregated[c(-2)]

paths_importance <- plot_misty_importance(MISTy_out = original_run_paths,
                                           pdf_out = "./revisions_misty/figures/importances_paths.pdf",
                                           predicted_features = predicted_paths,
                                           importance_cut = -10,
                                           midpoint = 0.5,
                                           make_clust = F,
                                           height_pdf = 4, 
                                           width_pdf = 3.5)

MISTy::plot_contrast_heatmap(misty.results = original_run,
                                            from.view = "intra",
                                            to.view = "para_path",
                                            cutoff = 0.5)


# Identifying importances of ligands ----------------------------------------------------------
predicted = original_run$improvements$target %>% unique()
original_run_ligands = original_run
original_run_ligands$importances.aggregated = original_run_ligands$importances.aggregated[c(-1,-3)]

# This is what is described in the manuscript ----------------------------------------------------------
ligands_importance = plot_misty_importance(MISTy_out = original_run_ligands,
                      pdf_out = "./revisions_misty/figures/importances_ligands.pdf",
                      predicted_features = predicted,
                      predictors_features = list(original_run$importances.aggregated$para_lig$Predictor),
                      importance_cut = 2,
                      make_clust = T,
                      midpoint = 0.5,
                      height_pdf = 4, 
                      width_pdf = 18)

ligands_importance_main = plot_misty_importance(MISTy_out = original_run_ligands,
                                           pdf_out = "./revisions_misty/figures/importances_ligands_main.pdf",
                                           predicted_features = predicted,
                                           predictors_features = list(original_run$importances.aggregated$para_lig$Predictor),
                                           importance_cut = 4.5,
                                           midpoint = 1,
                                           make_clust = T,
                                           height_pdf = 4, 
                                           width_pdf = 5)

ligands_importance = ligands_importance[[1]]

ligands_importance = ligands_importance %>% 
  mutate(Predictor = gsub("_","-",Predictor))

# Annotating each ligand and each receptor from the
# LR network coming from the subset of important ligands
# The output of this is the annotated information of each important ligand
ligands_importance_to_rec_anns = ligands_importance %>%
  group_by(Predicted) %>%
  dplyr::filter(Importance >= 2) %>%
  nest(lig_importances = c(Predictor,Importance)) %>% 
  mutate(regnetwork = map(lig_importances, function(x){ # Generate network including only important ligands per path
    
    receptor_network = lig_rec_red %>% 
      dplyr::filter(source_genesymbol %in% x$Predictor)
    
    return(receptor_network)
    
  })) %>% # Annotate Receptors
  mutate(receptor_ann = map(regnetwork, function(x){
    return(import_omnipath_annotations(proteins = unique(x$target_genesymbol)))
  })) %>% # Annotate Ligands
  mutate(ligand_ann = map(regnetwork, function(x){
    return(import_omnipath_annotations(proteins = unique(x$source_genesymbol)))
  })) %>% #Here you combine the annotations of receptors with the LR network
  dplyr::mutate(lig_interpretation = pmap(list(Predicted,regnetwork,receptor_ann), function(x,y,z){
    
    useful_anns = z %>%
      dplyr::filter(grepl(pattern = x, value, ignore.case = T))
    
    complete_ann = left_join(y,useful_anns, by = c("target_genesymbol" = "genesymbol"))
    
    complete_ann = complete_ann %>% 
      group_by(source_genesymbol) %>% 
      nest() %>%
      mutate(self_associated_receptors = map(data,na.omit)) %>%
      dplyr::select(-data)
    
    return(complete_ann)
    
  })) %>% #Here we finally see which ligands haven't been associated with a receptor annotated with the Predicted path
  dplyr::mutate(lig_importances = pmap(list(lig_importances,lig_interpretation), function(x,y){
    
    return(left_join(x,y,by = c("Predictor" = "source_genesymbol")))
    
  })) %>% left_join(view_contributions$data[[2]][,c("target","fraction")], 
                    by = c("Predicted"="target")) %>% 
  arrange(-fraction)


# Here I identify ligands that are part of the PROGENy response
# I assume that these ligands are by-product of pathway activation

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>% 
  mutate(interpretation = pmap(list(Predicted,lig_importances), function(x,y){
  
  path = gsub(pattern = "_",replacement = "-", x)
  
  pmat = getModel(organism = "Human", 1000) %>%
    rownames_to_column("gene") %>% pivot_longer(-gene,
                                                names_to = "pathway",
                                                values_to = "selfPROGENy_Score") %>%
    dplyr::filter(selfPROGENy_Score > 0)
  
  nov = left_join(y,pmat, by = c("Predictor"="gene"))
  
  return(nov)
}))

# Here additionally I will filter based on "important pathways for each marker"

para_path_contribution = view_contributions$data[[3]] %>%
  dplyr::select(target, fraction) %>%
  arrange(-fraction)

colnames(para_path_contribution) = c("target","para_path_fraction")

para_path_contribution = para_path_contribution %>%
  left_join(original_run$importances.aggregated$para_path %>% 
  pivot_longer(-Predictor,
               names_to = "Predicted",
               values_to = "Importance") %>%
  na.omit() %>%
  dplyr::filter(Importance > 1) %>%
  arrange(Predicted) %>%
  group_by(Predicted) %>% nest(path_importances = c(Predictor,Importance)), by = c("target" = "Predicted"))
  
# Final Annotation

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>%
  left_join(para_path_contribution, by = c("Predicted" = "target"))

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>% 
  dplyr::mutate(interpretation = pmap(list(Predicted, 
                                           interpretation, 
                                           path_importances, 
                                           lig_importances), 
                                        function(x,y,z,a){
    
    useful_p = c(x, z %>% pull(Predictor))
    
    interprt = y %>% group_by(Predictor) %>% nest()
    
    interprt = interprt %>% mutate(data = map(data, function(ann_lig,up){
      ann_lig_mod = ann_lig %>%
        mutate(pathway = ifelse(pathway %in% up,
                                pathway, NA)) %>%
        mutate(selfPROGENy_Score = ifelse(is.na(pathway),
                                NA, selfPROGENy_Score)) %>% 
        unique()
      
      return(ann_lig_mod)
    }, up = useful_p))
    
    progeny_info = interprt %>% 
      unnest() %>%
      ungroup() %>%
      dplyr::select(Predictor,pathway, selfPROGENy_Score) %>%
      group_by(Predictor) %>%
      nest(PROGENy_evidence = c(pathway, selfPROGENy_Score)) %>%
      dplyr::mutate(PROGENy_evidence = map(PROGENy_evidence, na.omit))
    
    final_interprt = left_join(a, progeny_info)
    
    return(final_interprt)
  }))

MISTy_results_annotation = ligands_importance_to_rec_anns %>% 
  dplyr::select(Predicted,
                interpretation,
                path_importances, 
                fraction,
                para_path_fraction)

saveRDS(MISTy_results_annotation, 
        file = "./revisions_misty/results/para_path_annotation.rds")

# These are the results are the ones presented in the manuscript -------------------
all_associations = MISTy_results_annotation %>%
  dplyr::select(Predicted, interpretation) %>%
  dplyr::mutate(interpretation_ns = map(interpretation, function(x){
    evidence = x %>%
      dplyr::select(-Importance) %>%
      dplyr::mutate(self_associated_receptors = map_dbl(self_associated_receptors, nrow),
                    PROGENy_evidence = map_dbl(PROGENy_evidence, nrow))
  })) %>% unnest(interpretation_ns)

all_associations %>% mutate(annotated = ifelse((self_associated_receptors == 0) &
                                               (PROGENy_evidence == 0), FALSE, TRUE)) %>%
  pull(annotated) %>% table()

### Plotting
# All useful ligands

important_ligands = MISTy_results_annotation %>% 
  dplyr::select(interpretation) %>% 
  unnest() %>%
  pull(Predictor) %>%
  unique()

important_ligands = gsub(pattern = "_",replacement = "-",important_ligands)

# Para matrices

para_ligs_2 = get_para_matrix(breast_A_1,para_assay = "SCT",l = 2,
                              para_features = important_ligands)

para_ligs_5 = get_para_matrix(breast_A_1,para_assay = "SCT",l = 5,
                              para_features = important_ligands)

para_ligs_10 = get_para_matrix(breast_A_1,para_assay = "SCT",l = 10,
                               para_features = important_ligands)


# Actually plotting
plot_para = function(pdf_file, width = 5, height = 3, para_mat, visium_slide){
  pdf(file = pdf_file,width = width,height = height)
  
  visium_slide[['para_mat']] = CreateAssayObject(data = para_mat)
  
  DefaultAssay(visium_slide) = "para_mat"
  
  for(feature in rownames(visium_slide)){
    plot(SpatialFeaturePlot(visium_slide,features = feature,stroke = 0))
  }
  
  dev.off()
}

plot_para(pdf_file = "./revisions_misty/paraligands_2.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_2)

plot_para(pdf_file = "./revisions_misty/paraligands_5.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_5)

plot_para(pdf_file = "./revisions_misty/paraligands_10.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_10)

plot_para(pdf_file = "./revisions_misty/paraligands_regularexpression.pdf", 
          visium_slide = breast_A_1, para_mat = as.matrix(breast_A_1@assays$SCT[important_ligands,]))

# Plotting receptors
receptors = MISTy_results_annotation %>% 
  dplyr::select(interpretation) %>% 
  unnest() %>% 
  dplyr::select(self_associated_receptors) %>% 
  unnest() %>%
  dplyr::filter(target_genesymbol %in% rownames(breast_A_1)) %>%
  pull(target_genesymbol) %>% 
  unique() %>% sort()

pdf(file = "./revisions_misty/receptors.pdf",width = 5,height = 3)

DefaultAssay(breast_A_1) = "SCT"

for(feature in receptors){
  plot(SpatialFeaturePlot(breast_A_1,features = feature,stroke = 0))
}

dev.off()

# Path

pdf(file = "./revisions_misty/pathways.pdf",width = 5,height = 3)
DefaultAssay(breast_A_1) = "progeny"
for(feature in rownames(breast_A_1)){
  plot(SpatialFeaturePlot(breast_A_1,features = feature,stroke = 0))
}
dev.off()

# Para path
para_ligs_2 = get_para_matrix(breast_A_1,para_assay = "progeny",l = 2)

para_ligs_5 = get_para_matrix(breast_A_1,para_assay = "progeny",l = 5)

para_ligs_10 = get_para_matrix(breast_A_1,para_assay = "progeny",l = 10)


plot_para(pdf_file = "./revisions_misty/paraoath_2.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_2)

plot_para(pdf_file = "./revisions_misty/parapath_5.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_5)

plot_para(pdf_file = "./revisions_misty/parapath_10.pdf", 
          visium_slide = breast_A_1, para_mat = para_ligs_10)


# Here auxiliary plots

progeny_mat = as.data.frame(t(as.matrix(breast_A_1@assays$progeny@data))) %>%
  rownames_to_column(var = "cell_id") %>% as_tibble() %>%
  pivot_longer(cell_id)

pdf("./revisions_misty/figures/intra_correlations.pdf", height = 4, width = 4)
print(ggplot(progeny_mat, aes(x = TNFa, y = NFkB)) +
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size =11)))

print(ggplot(progeny_mat, aes(x = MAPK, y = p53)) +
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size =11)))
dev.off()


DefaultAssay(breast_A_1) <- "progeny"

pdf(width = 4, height = 4.5, file = "./revisions_misty/figures/progeny_paths_main.pdf")

map(set_names(rownames(breast_A_1)), SpatialFeaturePlot, 
    object = breast_A_1, 
    stroke = 0,
    pt.size.factor = 2)

dev.off()

DefaultAssay(breast_A_1) <- "SCT"

pdf(width = 4, height = 4.5, file = "./revisions_misty/figures/genes_main.pdf")

map(c("STC1", "TNF", "TFF1", "EFNA1", "EDN1"), 
    SpatialFeaturePlot, 
    object = breast_A_1, 
    stroke = 0,
    pt.size.factor = 2)

dev.off()

