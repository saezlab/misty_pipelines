library(tidyverse)
library(MISTy)
library(cowplot)

bc.results <- collect_results(list.dirs("results/imc_bc/imc_bc_optim")[-1])
perm.results <- collect_results(seq(10) %>% map(~(list.dirs(paste0("results/imc_bc_perm", .x ,"/imc_bc_optim"))[-1])) %>% unlist)
meta <- read_delim("../Spatial/spatial_data/breastcancer/metadata.tsv", delim = "\t")

generate_plots <- function(global, meta){
  
  ymin <- min(global$value)
  ymax <- max(global$value)
  
  #global comparison
  global.plot <- ggplot(global, aes(x = set, y = value)) + 
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Global") +
    theme_classic()
  
  # grades
  
  grade1 <- global %>% filter(image %in% (meta %>% filter(Grade == 1) %>% pull(`Sample ID`)))
  grade1.plot <- ggplot(grade1, aes(x = set, y = value)) + 
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 1") +
    theme_classic()
  
  grade2 <- global %>% filter(image %in% (meta %>% filter(Grade == 2) %>% pull(`Sample ID`)))
  grade2.plot <- ggplot(grade2, aes(x = set, y = value)) + 
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 2") +
    theme_classic()
  
  grade3 <- global %>% filter(image %in% (meta %>% filter(Grade == 3) %>% pull(`Sample ID`)))
  grade3.plot <- ggplot(grade3, aes(x = set, y = value)) + 
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 3") +
    theme_classic()
  
  
  # tests
  
  
  list(global, grade1, grade2, grade3) %>% map_dbl(function(data){
    values <- data %>% group_by(set) %>% group_split() %>% map((. %>% pull(value)))
    wilcox.test(values[[1]], values[[2]], alternative = "greater")$p.value %>% log10
  })
  
  plot_grid(global.plot, grade1.plot, grade2.plot, grade3.plot)
  
}

global.perf <- bind_rows(
  tibble(set = "original", bc.results$improvements %>% filter(measure == "gain.R2") %>% 
           select(image, value)),
  tibble(set = "permutation", perm.results$improvements %>% filter(measure == "gain.R2") %>% 
           select(image, value))
) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))
generate_plots(global.perf, meta)


# plot coefficient fractions
global.intra <- bind_rows(
  tibble(set = "original", bc.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% ungroup() %>%
           filter(view == "intra") %>% select(image, value)),
  tibble(set = "permutation", perm.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>%  ungroup() %>%
           filter(view == "intra") %>% select(image, value))
) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))
generate_plots(global.intra, meta)


global.juxta <- bind_rows(
  tibble(set = "original", bc.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% ungroup() %>%
           filter(view == "juxta") %>% select(image, value)),
  tibble(set = "permutation", perm.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>%  ungroup() %>%
           filter(view == "juxta") %>% select(image, value))
) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))
generate_plots(global.juxta, meta)

global.para <- bind_rows(
  tibble(set = "original", bc.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>% ungroup() %>%
           filter(view == "para") %>% select(image, value)),
  tibble(set = "permutation", perm.results$contributions %>% filter(!str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
           group_by(image, target) %>% mutate(value = abs(value) / sum(abs(value))) %>%  ungroup() %>%
           filter(view == "para") %>% select(image, value))
) %>% mutate(image = str_extract(image, "[ABC][a-zA-Z0-9]+"))
generate_plots(global.para, meta)


frac.significant <- bind_rows(
  tibble(set = "original",
  bc.results$contributions %>% 
    filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
    group_by(view) %>% summarize(fraction = sum(value < 0.05)/n(), .groups = "drop")
  ),
  tibble(set = "permutation",
  perm.results$contributions %>% 
    filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>% 
    group_by(view) %>% summarize(fraction = sum(value < 0.05)/n(), .groups = "drop")
  )
)

ggplot(frac.significant, aes(x = view, y=fraction)) + 
  geom_point(aes(color = set), size = 2) + 
  ylim(0, 1) + theme_classic()
