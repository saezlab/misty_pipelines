# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Here we evaluate the footprint overlap between pathways

library(progeny)
library(tidyverse)

model <- progeny::getModel(organism = "Human", top = 1000) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  dplyr::filter(value != 0)

model_cor <- progeny::getModel(organism = "Human", top = 1000) %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column("path_a") %>%
  pivot_longer(-path_a) %>%
  dplyr::filter(path_a != name) %>%
  arrange(-abs(value))


model_cor$value %>% median()

pairwise <- progeny::getModel(organism = "Human", top = 1000)[, c("TNFa", "NFkB")] 
pairwise <- pairwise[which(rowSums(pairwise) != 0),]

tnfa_genes <- rownames(pairwise)[pairwise[, "TNFa"] != 0]
nfkb_genes <- rownames(pairwise)[pairwise[, "NFkB"] != 0]

length(intersect(nfkb_genes,tnfa_genes)) / length(union(nfkb_genes,tnfa_genes))
