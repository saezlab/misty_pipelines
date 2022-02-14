# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Interpretation of the original model

library(OmnipathR)
library(tidyverse)
library(progeny)
library(cowplot)

# Optimized aggregated performance --------------------------------------------
misty_results <- mistyR::collect_results(c("./results/misty_outs/A1/A1model_optim",
                                           "./results/misty_outs/A2/A2model_optim"))

para_path_contribution <- misty_results$contributions.stats %>%
  dplyr::select(target, view, mean) %>%
  dplyr::filter(view == "para_path") %>%
  dplyr::rename("para_path_fraction" = mean) %>%
  dplyr::select(-view)

para_lig_contribution <- misty_results$contributions.stats %>%
  dplyr::select(target, view, mean) %>%
  dplyr::filter(view == "para_lig") %>%
  dplyr::rename("fraction" = mean) %>%
  dplyr::select(-view)

# Importing ligand-receptor network --------------------------------------------
lig_rec <- import_intercell_network(
  interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
  transmitter_param = list(parent = 'ligand'),
  receiver_param = list(parent = 'receptor')
)

# Selecting important features --------------------------------------------
lig_rec_red <- lig_rec %>% dplyr::ungroup() %>%
  dplyr::select(target_genesymbol, source_genesymbol,
                is_directed, category_intercell_target)

# Ligands importances ------------------------------------------------------
ligands_importance <-  misty_results$importances.aggregated %>%
  dplyr::filter(view == "para_lig") %>%
  group_by(Predictor) %>%
  dplyr::mutate(max_importance = max(Importance)) %>%
  ungroup() %>%
  dplyr::filter(max_importance >= 2) 

ligands_importance <- ligands_importance %>% 
  mutate(Predictor = gsub("_","-",Predictor))

# Annotating each ligand and each receptor from the
# LR network coming from the subset of important ligands
# The output of this is the annotated information of each important ligand
ligands_importance_to_rec_anns <- ligands_importance %>%
  group_by(Target) %>%
  dplyr::select(-c("nsamples", "max_importance")) %>%
  dplyr::filter(Importance >= 2) %>%
  nest(lig_importances = c(Predictor,Importance)) %>% 
  mutate(regnetwork = map(lig_importances, function(x){ # Generate network including only important ligands per path
    
    receptor_network = lig_rec_red %>% 
      dplyr::filter(source_genesymbol %in% x$Predictor)
    
    return(receptor_network)
    
  })) %>% # Annotate Receptors
  mutate(receptor_ann = map(regnetwork, function(x){
    return(import_omnipath_annotations(proteins = unique(x$target_genesymbol)) %>%
             dplyr::filter(source != "PROGENy"))
  })) %>% # Annotate Ligands
  mutate(ligand_ann = map(regnetwork, function(x){
    return(import_omnipath_annotations(proteins = unique(x$source_genesymbol)) %>%
             dplyr::filter(source != "PROGENy"))
  })) %>% #Here you combine the annotations of receptors with the LR network
  dplyr::mutate(lig_interpretation = pmap(list(Target, regnetwork, receptor_ann), function(x,y,z){
    
    useful_anns = z %>%
      dplyr::filter(grepl(pattern = x, value, ignore.case = T))
    
    complete_ann = left_join(y, useful_anns, by = c("target_genesymbol" = "genesymbol"))
    
    complete_ann = complete_ann %>% 
      group_by(source_genesymbol) %>% 
      nest() %>%
      mutate(self_associated_receptors = map(data, na.omit)) %>%
      dplyr::select(-data)
    
    return(complete_ann)
    
  })) %>% #Here we finally see which ligands haven't been associated with a receptor annotated with the Predicted path
  dplyr::mutate(lig_importances = pmap(list(lig_importances, lig_interpretation), function(x,y){
    
    return(left_join(x,y,by = c("Predictor" = "source_genesymbol")))
    
  })) %>% left_join(para_lig_contribution, 
                    by = c("Target"="target")) %>% 
  arrange(-fraction)


# Here I identify ligands that are part of the PROGENy response
# I assume that these ligands are by-product of pathway activation

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>% 
  mutate(interpretation = pmap(list(Target, lig_importances), function(x,y){
    
    path = gsub(pattern = "_",replacement = "-", x)
    
    pmat = progeny::getModel(organism = "Human", 1000) %>%
      rownames_to_column("gene") %>% pivot_longer(-gene,
                                                  names_to = "pathway",
                                                  values_to = "selfPROGENy_Score") %>%
      dplyr::filter(selfPROGENy_Score > 0)
    
    nov = left_join(y,pmat, by = c("Predictor"="gene"))
    
    return(nov)
  }))

# Here additionally I will filter based on "important pathways for each marker"

para_path_contribution = para_path_contribution %>%
  left_join(misty_results$importances.aggregated %>%
              dplyr::filter(view == "para_path") %>%
              dplyr::select(-c("nsamples", "view")) %>%
              na.omit() %>%
              dplyr::filter(Importance > 1) %>%
              arrange(Target) %>%
              group_by(Target) %>% nest(path_importances = c(Predictor,Importance)), by = c("target" = "Target"))


# Final Annotation

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>%
  left_join(para_path_contribution, by = c("Target" = "target"))

ligands_importance_to_rec_anns = ligands_importance_to_rec_anns %>% 
  dplyr::mutate(interpretation = pmap(list(Target, 
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
  dplyr::select(Target,
                interpretation,
                path_importances, 
                fraction,
                para_path_fraction)



# These are the results are the ones presented in the manuscript -------------------
all_associations <- MISTy_results_annotation %>%
  dplyr::select(Target, interpretation) %>%
  dplyr::mutate(interpretation_ns = map(interpretation, function(x){
    evidence = x %>%
      dplyr::select(-Importance) %>%
      dplyr::mutate(self_associated_receptors = map(self_associated_receptors, nrow),
                    PROGENy_evidence = map_dbl(PROGENy_evidence, nrow)) %>%
      dplyr::mutate(self_associated_receptors = ifelse(is.null(self_associated_receptors), 0, self_associated_receptors))
  })) %>% unnest(interpretation_ns) %>%
  unnest(self_associated_receptors)

all_associations %>% mutate(annotated = ifelse((self_associated_receptors == 0) &
                                                 (PROGENy_evidence == 0), FALSE, TRUE)) %>%
  pull(annotated) %>% table()



