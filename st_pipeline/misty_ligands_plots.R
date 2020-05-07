# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de

#' This script generates the spatial transcriptomics figures
#' for main text after running misty_ligands.R

library(cowplot)
library(scales)
library(ggplot2)
source("slide_processing.R") #Check dependencies

# Define palettes
get_ligandf_colors = circlize::colorRamp2(seq(0,300,25),
                                          c("#003f5c",
                                            "#385b75",
                                            "#5f7a8f",
                                            "#8699aa",
                                            "#adbac6",
                                            "#d6dce2",
                                            "#ffffff",
                                            "#fde0e0",
                                            "#f9c2c1",
                                            "#f3a3a4",
                                            "#eb8387",
                                            "#e0636b",
                                            "#d43d51"))

get_pathf_colors = circlize::colorRamp2(seq(0,14,2),
                                        c("#003f5c",
                                          "#8699aa",
                                          "#d6dce2",
                                          "#ffffff",
                                          "#fde0e0",
                                          "#f3a3a4",
                                          "#eb8387",
                                          "#d43d51"))

SpatialPal = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))


breast_A_1 = readRDS(file = "breast_A_1.rds")
breast_A_2 = readRDS(file = "breast_A_2.rds")

# Generate output for paper
results_folder = c("mrun_breast1_optim",
                   "mrun_breast2_optim")

MISTy_out = MISTy_aggregator(results_folder = results_folder,
                             p.cutoff = 0.05)

# plots

# Supplemental 1-2. All pathways from both slides
DefaultAssay(breast_A_1) = "progeny"
allPROGENY_1 = SpatialPlot(breast_A_1, 
                           features = rownames(breast_A_1),
                           ncol = 4,
                           stroke = 0,
                           pt.size.factor = 2.5) +
  theme(legend.text = element_text(size=15),
        legend.title = element_text(size = 15))

pdf(file = "Sup_allPROGENY_1.pdf", width = 10,
    height = 10)

print(allPROGENY_1)

dev.off()


DefaultAssay(breast_A_2) = "progeny"
allPROGENY_2 = SpatialFeaturePlot(breast_A_2, 
                                  features = rownames(breast_A_2),
                                  ncol = 4,
                                  stroke = 0,
                                  pt.size.factor = 2.5) +
  theme(legend.text = element_text(size=15),
        legend.title = element_text(size = 15))

pdf(file = "Sup_allPROGENY_2.pdf", width = 10,
    height = 10)

print(allPROGENY_2)

dev.off()

# Supplemental 3. Improvement 
R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
  mutate(target = targets) %>%
  pivot_longer(cols = -target, 
               names_to = "name", 
               values_to = "value")

impr_plot = ggplot(R2_impr) +
  geom_point(aes(x = target, y = value * 100)) +
  theme_classic() +
  ylab("Change in variance explained") +
  xlab("Target") +
  #ylim(c(-5, 25)) +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))

plot(impr_plot)

# General Features Plot (Panel A)

path1 = SpatialPlot(breast_A_1, 
                    features = "nFeature_progeny",
                    ncol = 4,
                    stroke = 0,
                    pt.size.factor = 2.5) + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,
                                  size = 13.5),
        legend.text = element_text(size=13)) +
  labs(title="Number of pathways") +
  scale_fill_gradientn(colours = SpatialPal(length(seq(0,14,1))),
                       limits = c(0,14), 
                       breaks = c(0,5,10,14))
  

ligs1 = SpatialPlot(breast_A_1, 
                    features = "nFeature_ligands",
                    ncol = 4,
                    stroke = 0,
                    pt.size.factor = 2.5) + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,
                                  size = 13.5),
        legend.text = element_text(size=13)) +
  labs(title="Number of ligands") +
  scale_fill_gradientn(colours = SpatialPal(length(seq(100,275,1))),
                       limits = c(100,275), 
                       breaks = c(100,200,275))

title = ggdraw() + draw_label("Section 1", 
                              fontface='bold')

features_1 = plot_grid(path1, ligs1)

features_1 = plot_grid(title, features_1,
                       ncol=1, rel_heights=c(0.03, 1))

####

path2 = SpatialPlot(breast_A_2, 
                    features = "nFeature_progeny",
                    ncol = 4,
                    stroke = 0,
                    pt.size.factor = 2.5) + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,
                                  size = 13.5),
        legend.text = element_text(size=13)) +
  labs(title="Number of pathways") +
  scale_fill_gradientn(colours = SpatialPal(length(seq(0,14,1))),
                       limits = c(0,14), 
                       breaks = c(0,5,10,14))

ligs2 = SpatialPlot(breast_A_2, 
                    features = "nFeature_ligands",
                    ncol = 4,
                    stroke = 0,
                    pt.size.factor = 2.5) + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,
                                  size = 13),
        legend.text = element_text(size=13)) +
  labs(title="Number of ligands") +
  scale_fill_gradientn(colours = SpatialPal(length(seq(100,275,1))),
                       limits = c(100,275), 
                       breaks = c(100,200,275))

title = ggdraw() + draw_label("Section 2", 
                              fontface='bold')

features_2 = plot_grid(path2, ligs2)

features_2 = plot_grid(title, features_2,
                       ncol=1, rel_heights=c(0.03, 1))

#Contribution

coefs_plot = ggplot(MISTy_out$coefs) + 
  geom_col(aes(x=target, y=value, group=view, fill=view)) +
  xlab("Target") +
  ylab("Contribution") +
  theme_classic() +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5)) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("Intrinsic","Para pathway" ,"Para ligands")) +
  labs(fill = "View")


# First row

coefs_plot_mod = plot_grid(NULL, coefs_plot,
                          ncol=1, rel_heights=c(0.3, 1.2))

first_row = plot_grid(features_1,features_2,
                      coefs_plot_mod, nrow = 1,
                      rel_widths = c(.8,.8,1.3),
                      axis = "b")

pdf(file = "figure_st_up.pdf",
    width = 13,
    height = 4)


plot(first_row)


dev.off()


#Importance

#Intra pathways
importance_intra = tidyr::gather(MISTy_out$importance[[1]], "Predicted",
                                 "Importance", -Predictor)

as_tibble(importance_intra) %>% arrange(desc(Importance)) %>%
  slice(1:10)

importance_intra_plt = ggplot(importance_intra,
                              aes(x = Predictor, 
                                  y = Predicted, 
                                  fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size=10),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = scales::muted("blue"),
                       midpoint = 0.9) +
  xlab("Intrinsic pathways")

# Examples from the paper

DefaultAssay(breast_A_1) = "progeny"

intra_1 = SpatialFeaturePlot(object = breast_A_1,
                             features = "p53",
                             stroke = 0,
                             pt.size.factor = 2.) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(-4.5,4.5,0.25))),
                       limits = c(-4.5,4.5), 
                       breaks = c(-4,-2,0,2,4))

intra_2 = SpatialFeaturePlot(object = breast_A_1,
                             features = "MAPK",
                             stroke = 0,
                             pt.size.factor = 2.) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(-4.5,4.5,0.25))),
                       limits = c(-4.5,4.5), 
                       breaks = c(-4,-2,0,2,4))

intra_example = plot_grid(intra_1, intra_2,nrow = 1)
  
pdf(file = "intra_example.pdf",width = 5,height = 2.5)

plot(intra_example)

dev.off()


#Para_pathways
importance_para = tidyr::gather(MISTy_out$importance[[2]], "Predicted",
                                "Importance", -Predictor)

as_tibble(importance_para) %>% arrange(desc(Importance)) %>%
  slice(1:10)

importance_para_plt = ggplot(importance_para,
                             aes(x = Predictor, 
                                 y = Predicted, 
                                 fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size=10),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = scales::muted("blue"),
                       midpoint = 0.9) + 
  xlab("Para Pathways")

# Examples from the paper

DefaultAssay(breast_A_1) = "progeny"

para_1 = SpatialFeaturePlot(object = breast_A_1,
                             features = "JAK-STAT",
                             stroke = 0,
                             pt.size.factor = 2.) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(-4.5,4.5,0.25))),
                       limits = c(-6,6), 
                       breaks = c(-4,-2,0,2,4))

para_2 = SpatialFeaturePlot(object = breast_A_1,
                             features = "WNT",
                             stroke = 0,
                             pt.size.factor = 2.) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(-6,6,0.25))),
                       limits = c(-6,6), 
                       breaks = c(-4,-2,0,2,4))

parapath_example = plot_grid(para_1, para_2,nrow = 1)


pdf(file = "parapath_example.pdf",width = 5,height = 2.5)

plot(parapath_example)

dev.off()

#Ligands
importance_ligs = tidyr::gather(MISTy_out$importance[[3]], "Predicted",
                                "Importance", -Predictor)

informative_ligands = importance_ligs %>% group_by(Predictor) %>%
  summarise(mean_imp = mean(Importance)) %>%
  dplyr::filter(mean_imp > 2)

ligs_mat = MISTy_out$importance[[3]] %>% dplyr::filter(Predictor %in% informative_ligands$Predictor)
rownames(ligs_mat) = ligs_mat$Predictor
ligs_mat = as.matrix(ligs_mat[,-ncol(ligs_mat)])

ligs_clust =  hclust(d = t(dist(ligs_mat)))
path_clust =  hclust(d = t(dist(t(ligs_mat))))

importance_ligs_red = importance_ligs %>%
  dplyr::filter(Predictor %in% informative_ligands$Predictor) %>%
  dplyr::mutate(Predictor = factor(Predictor,
                                   levels = ligs_clust$labels[ligs_clust$order]),
                Predicted = factor(Predicted,
                                   levels = path_clust$labels[path_clust$order]))

importance_ligs_plt = ggplot(importance_ligs_red,
                             aes(x = Predictor, 
                                 y = Predicted, 
                                 fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        axis.text = element_text(size=10),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = scales::muted("blue"),
                       midpoint = 0.9) +
  xlab("Para Ligands")

# Examples from the paper

DefaultAssay(breast_A_1) = "SCT"

SEMA4B = SpatialFeaturePlot(object = breast_A_1,
                            features = c("SEMA4B"),
                            stroke = 0,
                            pt.size.factor = 2.5) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(0,2,0.25))),
                       limits = c(0,2), 
                       breaks = c(0,1,2))

DefaultAssay(breast_A_1) = "progeny"

Hypoxia = SpatialFeaturePlot(object = breast_A_1,
                             features = c("Hypoxia"),
                             stroke = 0,
                             pt.size.factor = 2.5) +
  scale_fill_gradientn(colours = SpatialPal(length(seq(-3,9,0.25))),
                       limits = c(-3,9), 
                       breaks = c(-3,0,3,6,9))

paralig_example = plot_grid(NULL,plot_grid(SEMA4B,Hypoxia, nrow = 1),NULL,
                            rel_widths = c(.1,1,.1),nrow = 1)

pdf(file = "paralig_example.pdf",width = 5,height = 2.5)

plot(paralig_example)

dev.off()

# Final plot

pdf(file = "figure_st_mid.pdf",
    width = 13,
    height = 4)

mid_panel = plot_grid(importance_intra_plt, importance_para_plt, 
                      importance_ligs_plt,
                      align = "h", nrow = 1,
                      rel_widths = c(0.6,0.6,1))

plot(mid_panel)


dev.off()





