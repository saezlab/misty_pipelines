# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Recovers ligands from Omnipath
library(tidyr)
library(tibble)
library(dplyr)
library(purrr)
library(OmnipathR)

InterCell_Annotations <- import_Omnipath_intercell() 
## We filter those proteins which are mainly annotated as receptor or ligand
Ligands_Receptors <- InterCell_Annotations %>%
  dplyr::filter(mainclass %in% c("receptor","ligand"))
## There are also some complexes. We are going to deal with them by including
## each of its individual proteins in our list
Ligand_Receptors_class <- character()
Ligand_Receptors_name <- character()
for (i in seq(nrow(Ligands_Receptors))){
  if (Ligands_Receptors$entity_type[i] == "complex"){
    Genescomplex <-unlist(strsplit(gsub("COMPLEX:", "", 
                                        Ligands_Receptors$genesymbol[i]),"_"))
    class <- rep(Ligands_Receptors$mainclass[i],length(Genescomplex))
    Ligand_Receptors_name <- c(Ligand_Receptors_name,Genescomplex)
    Ligand_Receptors_class <- c(Ligand_Receptors_class,class)
    
  } else {
    Ligand_Receptors_name <- 
      c(Ligand_Receptors_name, Ligands_Receptors$genesymbol[i]) 
    Ligand_Receptors_class <- 
      c(Ligand_Receptors_class, Ligands_Receptors$mainclass[i]) 
  }
}

Ligand_Receptors_df = data.frame(GeneSymbol = Ligand_Receptors_name, 
                                  Class = Ligand_Receptors_class, stringsAsFactors = FALSE) %>%
  dplyr::distinct()

AllLigands = dplyr::filter(Ligand_Receptors_df, Class == "ligand")

saveRDS(AllLigands,file = "opath_ligands.rds")
