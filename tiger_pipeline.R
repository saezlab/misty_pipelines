library(mistyR)
library(tidyverse)
library(future)
library(ggvoronoi)
library(cowplot)


# Run MISTy
data <- list.files("data/tiger", "\\.csv", full.names = TRUE)
plan(multisession)

l <- 15

data %>% walk(\(path){
  all <- read_csv(path)
  expr <- all %>% select(contains("gene"))
  pos <- all %>% select(x, y)
    
  expr.views <- create_initial_view(expr) %>% add_juxtaview(pos, l)
  run_misty(expr.views, bypass.intra = TRUE, 
            results.folder = paste0("results/tiger/expression/Tiger", 
                                    str_extract(path, "\\d")))
  
  ctype <- all %>% select(cell_id, cell_type) %>% mutate(value = 1) %>%
    group_by(cell_id) %>%
    pivot_wider(names_from = "cell_type", names_prefix = "ct", values_from = "value") %>%
    ungroup() %>%
    select(-cell_id) %>% replace(is.na(.),0)

  ctype.views <- create_initial_view(ctype) %>%
    add_juxtaview(pos, neighbor.thr = l, prefix = "juxta_")

  run_misty(ctype.views, bypass.intra = TRUE,
            results.folder = paste0("results/tiger/ctype/Tiger", 
                                    str_extract(path, "\\d")))
})

tiger2 <- read_csv(data[2]) %>% mutate(cell_type = as.factor(cell_type))

f2a <- ggplot(tiger2, aes(x,y, fill = cell_type )) + 
  geom_voronoi() + xlim(0,1000) + ylim(0,1000) + 
  theme_classic()

variances <- paste0("Tiger", seq_len(3)) %>% map_dfr(\(n){
  results  <- collect_results(paste0("results/tiger/ctype/", n))
  results$improvements %>% filter(measure == "gain.R2") %>%
    select(target, value) %>%
    mutate(sample = n)
})

f2b <- ggplot(variances, aes(x = sample, y = value, color = target)) + 
  geom_point() + 
  ylab("Variance explained") + 
  theme_classic()


tiger2.results <- collect_results("results/tiger/ctype/Tiger2/")
tiger2.results %>% plot_interaction_heatmap("juxta.15", cutoff = 0.5, trim = 5)
f2c1 <- last_plot() + ggtitle("Tiger 2 celltype juxtaview")

tiger3.results <- collect_results("results/tiger/ctype/Tiger3/")
tiger3.results %>% plot_interaction_heatmap("juxta.15", cutoff = 0.5, trim = 5)
f2c2 <- last_plot() + ggtitle("Tiger 3 celltype juxtaview")

tiger3.expr <- collect_results("results/tiger/expression/Tiger3")

tiger3.expr %>% plot_interaction_heatmap("juxta.15", clean = TRUE, cutoff = 2, trim = 1.8)
f2d <- last_plot() + ggtitle("Tiger 3 expression juxtaview")

plot_grid(plot_grid(f2a, f2b, f2c1, f2c2, labels = "AUTO"), f2d, rel_widths = c(4,1), labels = "AUTO")


