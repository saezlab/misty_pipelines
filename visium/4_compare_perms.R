# Copyright (c) [2022] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de

# Evaluate the performance of a MISTy model with a random layout

library(tidyverse)
library(mistyR)

# Optimized aggregated performance -------------------------------------------

original_run <- mistyR::collect_results(c("./results/misty_outs/A1/A1model_optim",
                                         "./results/misty_outs/A2/A2model_optim"))

# Permutation results ---------------------------------------------------------------

perm_files <- c(seq(5) %>% 
                 map(~ sprintf("./results/misty_outs/A1_rnd/A1model_rndm_%s_optim",.x)) 
               %>% unlist(),
               seq(5) %>% 
                 map(~ sprintf("./results/misty_outs/A2_rnd/A2model_rndm_%s_optim",.x)) %>% 
                 unlist())

perm_results <- collect_results(perm_files)

# 1. Compare R^2 of the multi-view models ---------------------------------------------

originalR2 = original_run$improvements %>% 
  dplyr::filter(measure == "multi.R2") %>%
  dplyr::select(target,value) %>%
  dplyr::mutate(layout = "original")

randomR2 = perm_results$improvements %>% 
  dplyr::filter(measure == "multi.R2") %>%
  dplyr::select(target,value) %>%
  dplyr::mutate(layout = "random")

R2_comp = ggplot(randomR2, aes(x = target, y = value, color = layout)) + 
  geom_boxplot() +
  geom_point(data = originalR2,aes(x = target, y = value, color = layout), 
             color = "navyblue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
        axis.text = element_text(size = 12)) +
  ylab("Model's R2")

pdf(file = "./results/figures/suppanel_a.pdf", height = 3, width = 3.7)

plot(R2_comp)

dev.off()

# Comparison of contributions and gain in importances ---------------------------
# Wilcoxon tests

generate_perm_plots <- function(global) {
  ymin <- min(global$value)
  ymax <- max(global$value)
  
  # global comparison
  global.plot <- ggplot(global, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Global") +
    theme_classic()
  
  # tests
  
  list(global) %>% map_dbl(function(data) {
    values <- data %>%
      group_by(set) %>%
      group_split() %>%
      map((. %>% pull(value)))
    wilcox.test(values[[1]], values[[2]], alternative = "greater")$p.value #%>% log10()
  }) %>% print
  
  return(global.plot)
  
}

generate_perm_plots_low <- function(global) {
  ymin <- min(global$value)
  ymax <- max(global$value)
  
  # global comparison
  global.plot <- ggplot(global, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Global") +
    theme_classic()
  
  # tests
  
  list(global) %>% map_dbl(function(data) {
    values <- data %>%
      group_by(set) %>%
      group_split() %>%
      map((. %>% pull(value)))
    wilcox.test(values[[1]], values[[2]], alternative = "lower")$p.value #%>% log10()
  }) %>% print
  
  return(global.plot)
  
}

# difference in gain.R2 -----------------------------------------------------------

global.perf <- bind_rows(
  tibble(set = "original", original_run$improvements %>% filter(measure == "gain.R2") %>%
           select(sample, value)),
  tibble(set = "permutation", perm_results$improvements %>% filter(measure == "gain.R2") %>%
           select(sample, value))
) %>% mutate(sample = str_extract(sample, "[ABC][a-zA-Z0-9]+"))

R2changepanel <- generate_perm_plots(global.perf)

R2changepanel <- R2changepanel + theme(axis.text = element_text(size = 12),
                                      axis.title.x = element_text(size = 15),
                                      axis.title.y = element_text(size = 15)) +
  ylab("Change in R2 (%)") +
  xlab("") + ggtitle("")

pdf(file = "./results/figures/suppanel_a_down.pdf", height = 3, width = 4)
plot(R2changepanel)  
dev.off()

# difference in contributions -------------------------------------------------

get_global_fractions <- function(bc.results, perm.results, from.view){
  bind_rows(
    tibble(set = "original", bc.results$contributions %>% 
             filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
             group_by(sample, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% 
             ungroup() %>%
             filter(view == from.view) %>% select(sample, value)),
    tibble(set = "permutation", perm.results$contributions %>% 
             filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
             group_by(sample, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% 
             ungroup() %>%
             filter(view == from.view) %>% select(sample, value))
  ) %>% mutate(sample = str_extract(sample, "[ABC][a-zA-Z0-9]+"))
}


global.intra  <- get_global_fractions(bc.results = original_run, perm.results = perm_results, "intra")
global.intra <- generate_perm_plots(global.intra)
global.intra <- global.intra + theme(axis.text = element_text(size = 12),
                                    axis.title.y = element_text(size = 15),
                                    legend.position = "none") +
  ylab("Intraview contribution") +
  xlab("") + 
  ggtitle("")

pdf(file = "./results/figures/suppanel_b_1.pdf", height = 3, width = 3)
plot(global.intra)  
dev.off()


global.paralig <- get_global_fractions(bc.results = original_run, perm.results = perm_results, "para_lig")
global.paralig <- generate_perm_plots(global.paralig)
global.paralig <- global.paralig + theme(axis.text = element_text(size = 12),
                                        axis.title.y = element_text(size = 15),
                                        legend.position = "none") +
  ylab("Ligand paraview \n contribution") +
  xlab("") + ggtitle("")

pdf(file = "./results/figures/suppanel_b_2.pdf", height = 3, width = 3)
plot(global.paralig)  
dev.off()

global.parapath = get_global_fractions(bc.results = original_run, perm.results = perm_results, "para_path")
global.parapath = generate_perm_plots(global.parapath)
global.parapath = global.parapath + theme(axis.text = element_text(size = 12),
                                          axis.title.y = element_text(size = 15),
                                          legend.position = "none") +
  ylab("Pathway paraview \n contribution") +
  xlab("") + ggtitle("")

pdf(file = "./results/figures/suppanel_b_3.pdf", height = 3, width = 3)
plot(global.parapath)  
dev.off()

# Panel --------
SF_viewchanges_plt <- cowplot::plot_grid(global.intra, 
                                         global.paralig, 
                                         global.parapath, nrow = 1, align = "hv")

pdf(file = "./results/figures/suppanel_b_all.pdf", height = 3, width = 9)
print(SF_viewchanges_plt)  
dev.off()









