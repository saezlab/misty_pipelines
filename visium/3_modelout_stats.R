# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Let's evaluate the model obtained in MISTy

library(mistyR)
library(tidyverse)

# Optimized aggregated performance --------------------------------------------
misty_results <- mistyR::collect_results(c("./results/misty_outs/A1/A1model_optim",
                                         "./results/misty_outs/A2/A2model_optim"))

# Improvement of R2 --------------------------------------------
measure <- "gain.R2"

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
                 axis.title.x = element_blank(),
                 axis.text = element_text(size = 11)) +
  ylab("Gain R2 (%)")

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
  theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     hjust = 1,
                                                     vjust = 0.5),
                 legend.position = "bottom",
                 axis.text = element_text(size = 11))

panel_b = cowplot::plot_grid(imp_results_plot, 
                             cont_results_plot, 
                             ncol = 1, 
                             align = "v", 
                             rel_heights = c(0.5, 1))

pdf(file = "./results/figures/panel_b.pdf", height = 6, width = 3)
plot(panel_b)
dev.off()

# Significance of change in R2 --------------------------------------------
# Adjusted p-values per image ---------------------------------------------
Rpvals_corrected <- misty_results$improvements %>%
  dplyr::filter(measure == "p.R2") %>% 
  group_by(sample) %>%
  dplyr::mutate(value = p.adjust(value)) %>% 
  ungroup()

RMSEpvals_corrected <- misty_results$improvements %>%
  dplyr::filter(measure == "p.RMSE") %>% 
  group_by(sample) %>%
  dplyr::mutate(value = p.adjust(value)) %>% 
  ungroup()

misty_results$improvements <- misty_results$improvements %>%
  dplyr::filter(measure != "p.RMSE",
                measure != "p.R2") %>%
  bind_rows(Rpvals_corrected,RMSEpvals_corrected)

#'Gets significant predictions from MISTy
bic <- misty_results$improvements %>% 
  dplyr::select(-sample) %>%
  group_by(target, measure) %>%
  summarise(mean_value = mean(value)) %>%
  pivot_wider(names_from = measure, values_from = mean_value) %>%
  dplyr::filter(p.R2 <= 0.1,
                  intra.R2 >= 0) %>% 
    select(target) %>% pull()

# Mean contribution of each view
misty_results$contributions.stats %>%
  group_by(view) %>%
  summarise(mean_contribution = mean(fraction))

# Quick summaries
improvement_summary <- misty_results$improvements.stats %>%
  group_by(measure) %>%
  arrange(-mean) %>%
  nest()

view_contributions <- misty_results$contributions.stats %>%
  group_by(view) %>%
  arrange(-mean) %>%
  nest()

# Text observations
view_contributions %>%
  dplyr::filter(view == "para_path") %>%
  unnest()
para_path

# Plot importances:

# Intra

intra_plt <- misty_results$importances.aggregated %>%
  dplyr::filter(view == "intra") %>%
  dplyr::select(-nsamples) %>%
  ggplot(aes(x = Predictor, 
             y = Target,
             fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size=11),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = "#8DA0CB",
                       midpoint = 0) +
  xlab("Predictors") + ylab("Targets") +
  coord_equal()

# Para path

para_plt <- misty_results$importances.aggregated %>%
  dplyr::filter(view == "para_path") %>%
  dplyr::select(-nsamples) %>%
  ggplot(aes(x = Predictor, 
             y = Target,
             fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size=11),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = "#8DA0CB",
                       midpoint = 0) +
  xlab("Predictors") + ylab("Targets") +
  coord_equal()


pdf(file = "./results/figures/panel_c.pdf", height = 4, width = 4.5)
plot(intra_plt)
dev.off()

pdf(file = "./results/figures/panel_e.pdf", height = 4, width = 4.5)
plot(para_plt)
dev.off()

# Ligand plot

# First identify ligands >=2 importance

filtered_ligands <-  misty_results$importances.aggregated %>%
  dplyr::filter(view == "para_lig") %>%
  group_by(Predictor) %>%
  dplyr::mutate(max_importance = max(Importance)) %>%
  ungroup() %>%
  dplyr::filter(max_importance >= 2.5) 
  
filtered_ligands %>% pull(Predictor) %>% unique() %>% length()

filtered_ligands_mat <- filtered_ligands %>%
  dplyr::select(Predictor, Target, Importance) %>%
  pivot_wider(names_from = Target, values_from = Importance) %>%
  column_to_rownames("Predictor")

predictor_order <- hclust(dist(filtered_ligands_mat))
predictor_order <- predictor_order$labels[predictor_order$order]

target_order <- hclust(dist(filtered_ligands_mat %>% t()))
target_order <- target_order$labels[target_order$order]

para_lig_plt <- filtered_ligands %>%
  dplyr::mutate(Predictor = factor(Predictor, levels = predictor_order),
                Target = factor(Target, levels = target_order)) %>%
  ggplot(aes(x = Predictor, 
            y = Target,
            fill = Importance)) + geom_tile() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size=11),
        legend.key.size = unit(.6, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "white", 
                       mid = "white", 
                       high = "#8DA0CB",
                       midpoint = 0) +
  xlab("Predictors") + ylab("Targets") +
  coord_equal()

pdf(file = "./results/figures/suppanel_d.pdf", height = 4, width = 16)
plot(para_lig_plt)
dev.off()
