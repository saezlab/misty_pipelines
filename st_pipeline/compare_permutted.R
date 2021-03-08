# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de

#' Compares permuted analysis with original one:
#' -Shows the increment of relevance of the para-view in original data
#' -Shows that the change in R2 is not random
#' 

library(OmnipathR)
library(tidyverse)
library(Seurat)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")


# Optimized aggregated performance -------------------------------------------

original_run = MISTy::collect_results(c("./revisions_misty/results/A1model_optim",
                                        "./revisions_misty/results/A2model_optim"))


original_run %>% plot_improvement_stats() %>% plot_view_contributions()


# Permutation results ---------------------------------------------------------------

perm_files = c(seq(5) %>% 
                 map(~ sprintf("./revisions_misty/results/A1model_rndm_%s_optim",.x)) 
               %>% unlist(),
               seq(5) %>% 
                 map(~ sprintf("./revisions_misty/results/A2model_rndm_%s_optim",.x)) %>% 
                 unlist())

perm_results = collect_results(perm_files)

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
  geom_point(data = originalR2,aes(x = target, y = value, color = layout), color = "lightblue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
        axis.text = element_text(size = 11)) +
  ylab("R^2")

pdf(file = "./revisions_misty/figures/SF_R2panel.pdf", height = 5, width = 6)

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

# difference in gain.R2 -----------------------------------------------------------

global.perf <- bind_rows(
  tibble(set = "original", original_run$improvements %>% filter(measure == "gain.R2") %>%
           select(image, value)),
  tibble(set = "permutation", perm_results$improvements %>% filter(measure == "gain.R2") %>%
           select(image, value))
) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))

R2changepanel = generate_perm_plots(global.perf)

R2changepanel = R2changepanel + theme(axis.text = element_text(size = 11)) +
  ylab("Change in R^2 (in percentage)") +
  xlab("") + ggtitle("")

pdf(file = "./revisions_misty/figures/SF_R2changepanel.pdf", height = 5, width = 6)
plot(R2changepanel)  
dev.off()

# difference in contributions -------------------------------------------------

get_global_fractions <- function(bc.results, perm.results, from.view){
  bind_rows(
    tibble(set = "original", bc.results$contributions %>% 
             filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
             group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% 
             ungroup() %>%
             filter(view == from.view) %>% select(image, value)),
    tibble(set = "permutation", perm.results$contributions %>% 
             filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
             group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% 
             ungroup() %>%
             filter(view == from.view) %>% select(image, value))
  ) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))
}


global.intra  = get_global_fractions(bc.results = original_run, perm.results = perm_results, "intra")
global.intra = generate_perm_plots(global.intra)
global.intra = global.intra + theme(axis.text = element_text(size = 11)) +
  ylab("Intraview contribution") +
  xlab("") + ggtitle("")

pdf(file = "./revisions_misty/figures/SF_globalintrapanel.pdf", height = 5, width = 6)
plot(global.intra)  
dev.off()


global.paralig = get_global_fractions(bc.results = original_run, perm.results = perm_results, "para_lig")
global.paralig = generate_perm_plots(global.paralig)
global.paralig = global.paralig + theme(axis.text = element_text(size = 11)) +
  ylab("Ligand paraview contribution") +
  xlab("") + ggtitle("")

pdf(file = "./revisions_misty/figures/SF_globalparaligpanel.pdf", height = 5, width = 6)
plot(global.paralig)  
dev.off()

global.parapath = get_global_fractions(bc.results = original_run, perm.results = perm_results, "para_path")
global.parapath = generate_perm_plots(global.parapath)
global.parapath = global.parapath + theme(axis.text = element_text(size = 11)) +
  ylab("Pathway paraview contribution") +
  xlab("") + ggtitle("")

pdf(file = "./revisions_misty/figures/SF_globalparapathpanel.pdf", height = 5, width = 6)
plot(global.parapath)  
dev.off()


frac.significant <- bind_rows(
  tibble(
    set = "original",
    original_run$contributions %>%
      filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
      group_by(view) %>% summarize(fraction = sum(value < 0.05) / n(), .groups = "drop")
  ),
  tibble(
    set = "permutation",
    perm_results$contributions %>%
      filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
      group_by(view) %>% summarize(fraction = sum(value < 0.05) / n(), .groups = "drop")
  )
)

ggplot(frac.significant, aes(x = view, y = fraction)) +
  geom_point(aes(color = set), size = 2) +
  ylim(0, 1) +
  theme_classic()

