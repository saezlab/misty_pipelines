library(mistyR)
library(tidyverse)
library(future)
library(ggvoronoi)
library(cowplot)


# Run MISTy
data <- list.files("data/insilico/structure", "tissue.*\\.csv", full.names = TRUE)
plan(multisession)

l <- 15

data %>% walk(\(path){
  all <- read_csv(path)
  expr <- all %>% select(contains("gene"))
  pos <- all %>% select(x, y)
    
  expr.views <- create_initial_view(expr) %>% add_juxtaview(pos, l)
  run_misty(expr.views, bypass.intra = TRUE, 
            results.folder = paste0("results/structure/expression/tissue", 
                                    str_extract(path, "\\d")))
  
  ctype <- all %>% select(cell_id, cell_type) %>% mutate(value = 1) %>%
    group_by(cell_id) %>%
    pivot_wider(names_from = "cell_type", names_prefix = "ct", values_from = "value") %>%
    ungroup() %>%
    select(-cell_id) %>% replace(is.na(.),0)

  ctype.views <- create_initial_view(ctype) %>%
    add_juxtaview(pos, neighbor.thr = l, prefix = "juxta_")

  run_misty(ctype.views, bypass.intra = TRUE,
            results.folder = paste0("results/structure/ctype/tissue", 
                                    str_extract(path, "\\d")))
})

tissues <- data %>% map_dfr(~ 
  read_csv(.x) %>% 
    mutate(cell_type = as.factor(cell_type), 
           tissue=paste("Tissue", str_extract(.x, "\\d")))
)

f2a <- ggplot(tissues, aes(x,y, fill = cell_type )) + 
geom_voronoi() + xlim(0,1000) + ylim(0,1000) + 
coord_fixed() +
facet_wrap(vars(tissue)) +
theme_classic()

variances <- paste0("tissue", seq_len(3)) %>% map_dfr(\(n){
  results  <- collect_results(paste0("results/structure/ctype/", n))
  results$improvements %>% filter(measure == "gain.R2") %>%
    select(target, value) %>%
    mutate(sample = n)
})

f2b <- ggplot(variances, aes(x = sample, y = value, color = target)) + 
  geom_point() + 
  ylab("Variance explained") + 
  theme_classic()


tissue2.results <- collect_results("results/structure/ctype/tissue2/")
tissue2.results %>% plot_interaction_heatmap("juxta.15", cutoff = 0.5, trim = 5)
f2c1 <- last_plot() + ggtitle("Tissue 2 celltype juxtaview")

tissue3.results <- collect_results("results/structure/ctype/tissue3/")
tissue3.results %>% plot_interaction_heatmap("juxta.15", cutoff = 0.5, trim = 5)
f2c2 <- last_plot() + ggtitle("Tissue 3 celltype juxtaview")

tissue3.results.expr <- collect_results("results/structure/expression/tissue3")

tissue3.results.expr %>% plot_interaction_heatmap("juxta.15", clean = TRUE, 
                                                  cutoff = 2, trim = 4)
f2d <- last_plot() + ggtitle("Tissue 3 expression juxtaview")

# Predictor topmarker importances

tissue3 <- read.csv(data[3])

ct0.expr <- tissue3 %>% filter(cell_type == 0) %>% select(-seq_len(5)) %>% colMeans
non.ct0.expr <- tissue3 %>% filter(cell_type != 0) %>% select(-seq_len(5)) %>% colMeans
ct0.markers <- ct0.expr - non.ct0.expr

ct1.expr <- tissue3 %>% filter(cell_type == 1) %>% select(-seq_len(5)) %>% colMeans
non.ct1.expr <- tissue3 %>% filter(cell_type != 1) %>% select(-seq_len(5)) %>% colMeans
ct1.markers <- ct1.expr - non.ct1.expr


ct2.expr <- tissue3 %>% filter(cell_type == 2) %>% select(-seq_len(5)) %>% colMeans
non.ct2.expr <- tissue3 %>% filter(cell_type != 2) %>% select(-seq_len(5)) %>% colMeans
ct2.markers <- ct2.expr - non.ct2.expr

ct3.expr <- tissue3 %>% filter(cell_type == 3) %>% select(-seq_len(5)) %>% colMeans
non.ct3.expr <- tissue3 %>% filter(cell_type != 3) %>% select(-seq_len(5)) %>% colMeans
ct3.markers <- ct3.expr - non.ct3.expr

topn <- 10

top0 <- names(sort(abs(ct0.markers), decreasing = TRUE))[1:topn]
top1 <- names(sort(abs(ct1.markers), decreasing = TRUE))[1:topn]
top2 <- names(sort(abs(ct2.markers), decreasing = TRUE))[1:topn]
top3 <- names(sort(abs(ct3.markers), decreasing = TRUE))[1:topn]

markers <- data.frame(cell_type = paste0("ct", rep(seq(4)-1, 10)) %>% sort(), marker = c(top0, top1, top2, top3))


cts <- c("ct0", "ct2")

importances <- expand_grid(cts, markers %>% pull(cell_type) %>% unique()) %>% pmap_dfr(~
data.frame(Target = ..1, Predictor = ..2, Importance = tissue3.results.expr$importances.aggregated %>%
  filter(
    Target %in% (markers %>% filter(cell_type == ..1) %>% pull(marker) %>% 
                   intersect(tissue3.results.expr$improvements.stats %>% 
                               filter(measure == "gain.R2", mean >= 4) %>% 
                               pull(target))),
    Predictor %in% (markers %>% filter(cell_type == ..2) %>% pull(marker))
  ) %>%
  pull(Importance)) %>% drop_na())


f2e <- ggplot(importances, aes(x = Predictor, y = Importance, fill = Predictor)) + 
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_hline(yintercept = 0, col = "gray30", linetype = "dotted") +
  facet_wrap(vars(Target), ncol = 1) + theme_classic() +
  theme(legend.position = "none")


pdf("plots/structure/Figure2.pdf", width = 13.2, height = 16.6)
plot_grid(
plot_grid(
  f2a, f2b, nrow = 1, rel_widths = c(1.75,1), labels = c("A","B")
),
plot_grid(plot_grid(f2c1, f2c2, ncol = 1), f2d, f2e, nrow = 1, rel_widths = c(0.8,1,0.8), labels = c("C", "D", "E")),
ncol = 1,
rel_heights = c(1, 1.5)
)
dev.off()

