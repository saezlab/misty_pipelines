library(tidyverse)
library(magrittr)
library(mistyR)
library(cowplot)
library(factoextra)

future::plan(future::multisession)

bc.results <- collect_results(list.dirs("results/imc_small_ridge_perm0/imc_bc_optim_zoi")[-1])
meta <- read_delim("data/imc_small_breastcancer/metadata.tsv", delim = "\t") %>%
  filter(`Sample ID` %in% list.dirs("results/imc_small_ridge_perm0/imc_bc_optim_zoi",
    full.names = FALSE
  ))


bc.results %<>% map(~ .x %>% mutate(across(
  matches("([tT]arget|Predictor)"),
  \(x) str_replace_all(x, "Cytokeratin", "CK") %>%
    str_replace("CarbonicAnhydrase", "CA") %>%
    str_replace("CleavedCaspase", "CC") %>%
    str_replace("Histone", "H")
)))


bc.results %>% plot_improvement_stats()

f3a <- last_plot()

bc.results %>% plot_improvement_stats("intra.R2")

sf3a1 <- last_plot()

bc.results %>% plot_improvement_stats("multi.R2")

sf3a2 <- last_plot()


bc.results %>% plot_view_contributions(trim = 1)

f3b <- last_plot()

f3b.extra <- bc.results$contributions.stats %>%
  filter(target %in%
    (bc.results$improvements.stats %>%
      filter(measure == "gain.R2") %>%
      filter(mean >= 1) %>% pull(target))) %>%
  group_by(view) %>%
  summarise(mfrac = mean(fraction)) %>%
  ggplot(aes(x = "    Average", y = mfrac, group = view, fill = view)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  xlab("") +
  ylab("Contribution") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

g1.folders <- paste0(
  "results/imc_small_ridge_perm0/imc_bc_optim_zoi/",
  meta %>% filter(Grade == 1) %>% pull(`Sample ID`)
)

g3.folders <- paste0(
  "results/imc_small_perm0/imc_bc_optim/",
  meta %>% filter(Grade == 3) %>% pull(`Sample ID`)
)


grade1.results <- collect_results(g1.folders)
grade1.results %<>% map(~ .x %>% mutate(across(
  matches("([tT]arget|Predictor)"),
  \(x) str_replace_all(x, "Cytokeratin", "CK") %>%
    str_replace("CarbonicAnhydrase", "CA") %>%
    str_replace("CleavedCaspase", "CC") %>%
    str_replace("Histone", "H")
)))

grade3.results <- collect_results(g3.folders)
grade3.results %<>% map(~ .x %>% mutate(across(
  matches("([tT]arget|Predictor)"),
  \(x) str_replace_all(x, "Cytokeratin", "CK") %>%
    str_replace("CarbonicAnhydrase", "CA") %>%
    str_replace("CleavedCaspase", "CC") %>%
    str_replace("Histone", "H")
)))

grade1.results %>% plot_contrast_heatmap("intra", "para", 0.5, trim = 1)
f4b1 <- last_plot() + ggtitle("Grade 1 paraview - intraview")

grade3.results %>% plot_contrast_heatmap("intra", "para", 0.5, trim = 1)
f4b2 <- last_plot() + ggtitle("Grade 3 paraview - intraview")

grade3.results %>% 
  plot_contrast_results(grade1.results, views = "intra", 
                        cutoff.from = 0.5, cutoff.to = 0.5, trim = 1)
f4d1 <- last_plot() + ggtitle("Intraview Grade 1 - Grade 3")

grade3.results %>% 
  plot_contrast_results(grade1.results, views = "para", 
                        cutoff.from = 0.5, cutoff.to = 0.5, trim = 1)
f4d2 <- last_plot() + ggtitle("Paraview Grade 1 - Grade 3")


# Reduction in important interactions

imp.inter <- rbind(
  grade1.results$importances.aggregated %>%
    filter(Importance >= 0.5) %>%
    group_by(view) %>%
    summarize(value = n()) %>%
    mutate(grade = "grade 1"),
  grade3.results$importances.aggregated %>%
    filter(Importance >= 0.5) %>%
    group_by(view) %>%
    summarize(value = n()) %>%
    mutate(grade = "grade 3")
)

f4c <- ggplot(imp.inter, aes(x = grade, y = value, group = view, color = view)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  xlab("Grade") +
  ylab("Important interactions") +
  ylim(40, 100) +
  theme_classic()

# Permutation analysis

perm.results <- collect_results(seq(10) %>%
  map(~ (list.dirs(paste0("results/imc_small_ridge_perm", .x, "/imc_bc_optim_zoi"))[-1])) %>%
  unlist())


perm.results %<>% map(~ .x %>% mutate(across(
  matches("([tT]arget|Predictor)"),
  \(x) str_replace_all(x, "Cytokeratin", "CK") %>%
    str_replace("CarbonicAnhydrase", "CA") %>%
    str_replace("CleavedCaspase", "CC") %>%
    str_replace("Histone", "H")
)))


perm.results %>% plot_improvement_stats()
sf3b1 <- last_plot()

perm.results %>% plot_view_contributions()
sf3b2 <- last_plot()


get_global_fractions <- function(bc.results, perm.results, from.view) {
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


# difference in gain

global.perf <- bind_rows(
  tibble(set = "original", bc.results$improvements %>% filter(measure == "gain.R2") %>%
    select(sample, value)),
  tibble(set = "permutation", perm.results$improvements %>% filter(measure == "gain.R2") %>%
    select(sample, value))
) %>% mutate(sample = str_extract(sample, "[ABC][a-zA-Z0-9]+"))


ymin <- min(global.perf$value)
ymax <- max(global.perf$value)

# global comparison
f3c <- ggplot(global.perf, aes(x = set, y = value)) +
  geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set2") +
  ylim(ymin, ymax) +
  ggtitle("Gain R2") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())


# plot coefficient fractions

global.intra <- get_global_fractions(bc.results, perm.results, "intra")

ymin <- min(global.intra$value)
ymax <- max(global.intra$value)

f3d1 <- ggplot(global.intra, aes(x = set, y = value)) +
  geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set2") +
  ylim(ymin, ymax) +
  ggtitle("Intraview contribution") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())


global.juxta <- get_global_fractions(bc.results, perm.results, "juxta")

ymin <- min(global.juxta$value)
ymax <- max(global.juxta$value)

f3d2 <- ggplot(global.juxta, aes(x = set, y = value)) +
  geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set2") +
  ylim(ymin, ymax) +
  ggtitle("Juxtaview contribution") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())

global.para <- get_global_fractions(bc.results, perm.results, "para")

ymin <- min(global.para$value)
ymax <- max(global.para$value)

f3d3 <- ggplot(global.para, aes(x = set, y = value)) +
  geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set2") +
  ylim(ymin, ymax) +
  ggtitle("Paraview contribution") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())


# Singnature analysis

signature <- extract_signature(bc.results, trim = 1)

signature.pca <- prcomp(signature %>% select(-sample))

meta.pca <- left_join(as_tibble(signature.pca$x) %>%
  mutate(sample = signature %>% pull(sample) %>% str_extract("[ABC]y[0-9]+x[0-9]+")),
meta %>% filter(`Sample ID` %in% (signature %>% pull(sample) %>% str_extract("[ABC]y[0-9]+x[0-9]+"))),
by = c("sample" = "Sample ID")
) %>%
  mutate(Grade = as.factor(Grade))

f3e1 <- ggplot(meta.pca %>% filter(!is.na(Grade)), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Grade), size = 3) +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

f3e2 <- fviz_pca_var(signature.pca,
  col.var = "cos2", repel = TRUE, select.var = list(cos2 = 10),
  gradient.cols = c("#666666", "#377EB8", "#E41A1C")
) + theme_classic()

# generate importance signature

imp.signature <- extract_signature(bc.results, type = "importance", trim = 1)

impsig.pca.raw <- prcomp(imp.signature %>% select(-sample))

impsig.pca <- data.frame(sample = bc.results$improvements %>% pull(sample) %>% unique() %>%
  str_extract("[ABC][a-zA-Z0-9]+"), impsig.pca.raw$x)

impmeta.pca <- left_join(meta %>% filter(`Sample ID` %in% (impsig.pca %>% pull(sample))),
  impsig.pca,
  by = c("Sample ID" = "sample")
) %>% mutate(Grade = as.factor(Grade))

f4a1 <- ggplot(impmeta.pca %>% filter(!is.na(Grade)), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Grade), size = 3) +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

f4a2 <- fviz_pca_var(impsig.pca.raw,
  col.var = "cos2", select.var = list(cos2 = 10),
  gradient.cols = c("#666666", "#377EB8", "#E41A1C"), repel = TRUE
) + theme_classic()


# raw expression

data <- list.dirs("data/imc_small_breastcancer/", recursive = FALSE)

mean.expr <- data %>% map_dfr(function(d) {
  expr <- read_delim(paste0(d, "/expressions.txt"), delim = " ", col_types = cols())
  colnames(expr) <- make.names(colnames(expr))
  c(image = str_extract(d, "[ABC][a-zA-Z0-9]+"), colMeans(expr))
})

colnames(mean.expr) %<>% str_replace_all("Cytokeratin", "CK") %>%
  str_replace("CarbonicAnhydrase", "CA") %>%
  str_replace("CleavedCaspase", "CC") %>%
  str_replace("Histone", "H")

mean.expr.pca <- mean.expr %>%
  select(-image) %>%
  mutate_all(as.numeric) %>%
  prcomp(scale. = TRUE)


mean.expr.join <- left_join(cbind(mean.expr, mean.expr.pca$x), meta, by = c("image" = "Sample ID")) %>%
  mutate(Grade = as.factor(Grade)) %>%
  filter(!is.na(Grade))

sf3c1 <- ggplot(mean.expr.join, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Grade), size = 3) +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

sf3c2 <- fviz_pca_var(mean.expr.pca,
  col.var = "cos2", repel = TRUE, col.circle = NA,
  gradient.cols = c("#666666", "#377EB8", "#E41A1C")
) + theme_classic()


# Figure 3
pdf("plots/imc_small/Figure3.pdf", width = 13.2, height = 16.6)
plot_grid(f3a,
  plot_grid(f3b, f3b.extra, rel_widths = c(4, 1), axis = "b", align = "v"),
  f3e1, f3e2, f3c,
  f3d1, f3d2, f3d3,
  ncol = 2, byrow = FALSE, rel_widths = c(2, 1),
  labels = c("A", "C", "B", "D", "E", "", "", "")
)
dev.off()

# Figure 4
pdf("plots/imc_small/Figure4.pdf", width = 13.2, height = 16.6)
plot_grid(plot_grid(f4a1, f4a2, f4b1, f4b2,
  ncol = 2, byrow = TRUE,
  scale = c(0.8, 1, 1, 1), rel_heights = c(1.5, 1),
  labels = c("A", "", "B", "")
),
plot_grid(f4c, f4d1, f4d2,
  nrow = 1,
  labels = c("C", "D", "")
),
ncol = 1, rel_heights = c(2, 1)
)
dev.off()

# Supplementary Figure 3
pdf("plots/imc_small/SupplementaryFigure3.pdf", width = 13.2, height = 16.6)
plot_grid(sf3a1, sf3a2, sf3b1, sf3b2, sf3c1, sf3c2,
          ncol = 2, byrow = TRUE,
          labels = c("A","", "B", "", "C", ""))
dev.off()
