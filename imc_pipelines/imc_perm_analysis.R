library(tidyverse)
library(mistyR)
library(cowplot)
library(factoextra)

future::plan(future::multisession)


bc.results <- collect_results(list.dirs("results/imc_small_ridge_perm0/imc_bc_optim_zoi")[-1])
meta <- read_delim("data/imc_small_breastcancer/metadata.tsv", delim = "\t") %>% 
  filter(`Sample ID` %in% list.dirs("results/imc_small_ridge_perm0/imc_bc_optim_zoi",  
                                    full.names = FALSE))

# Plot results all samples

# Figure S3
bc.results %>% plot_improvement_stats()

# Figure 3
bc.results %>%
  plot_view_contributions(trim = 1) %>%
  plot_interaction_heatmap("intra", 0.5) %>%
  plot_interaction_heatmap("juxta", 0.5, trim = 1) %>%
  plot_interaction_heatmap("para", 0.5, trim = 1) %>%
  plot_interaction_communities("intra", 0.5) %>%
  plot_interaction_communities("juxta", 0.5) %>%
  plot_interaction_communities("para", 0.5)


# Plots per grade
plot_grade_collection <- function(grade.results, grade){
  pdf(paste0("plots/imc_small/grade_", grade ,".pdf"), width = 6.5, height = 5)
  grade.results %>%
    plot_improvement_stats() %>%
    plot_view_contributions(trim = 1) %>%
    plot_interaction_heatmap("intra", 0.5) %>%
    plot_interaction_heatmap("juxta", 0.5, trim = 1) %>%
    plot_interaction_heatmap("para", 0.5, trim = 1) %>%
    plot_interaction_communities("intra", 0.5) %>%
    plot_interaction_communities("juxta", 0.5) %>%
    plot_interaction_communities("para", 0.5) %>%
    plot_contrast_heatmap("intra", "juxta", 0.5, trim = 1) %>%
    plot_contrast_heatmap("intra", "para", 0.5, trim = 1)
  
  dev.off()
}

ggplot(bc.results$contributions.stats %>% 
         group_by(view) %>% 
         summarise(mfrac = mean(fraction)), 
       aes(x = "", y = mfrac, group = view, fill = view)) + 
  geom_col() + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Average") +
  ylab("Contribution") +
  theme_classic()
  

g1.folders <- paste0("results/imc_small_ridge_perm0/imc_bc_optim_zoi/", 
                     meta %>% filter(Grade == 1) %>% pull(`Sample ID`))
g2.folders <- paste0("results/imc_small_perm0/imc_bc_optim/", 
                     meta %>% filter(Grade == 2) %>% pull(`Sample ID`))
g3.folders <- paste0("results/imc_small_perm0/imc_bc_optim/", 
                     meta %>% filter(Grade == 3) %>% pull(`Sample ID`))


grade1.results <- collect_results(g1.folders)
grade2.results <- collect_results(g2.folders)
grade3.results <- collect_results(g3.folders)

plot_grade_collection(grade1.results, 1)
plot_grade_collection(grade2.results, 2)
plot_grade_collection(grade3.results, 3)

grade1.results %>% plot_contrast_heatmap("intra", "para", 0.5, trim = 1)
grade2.results %>% plot_contrast_heatmap("intra", "para", 0.5, trim = 1)
grade3.results %>% plot_contrast_heatmap("intra", "para", 0.5, trim = 1)

grade2.results %>% plot_contrast_results(grade1.results, cutoff.from = 1, cutoff.to = 1)

pdf(file = "plots/imc_small/small_contrast_g3g1.pdf", width = 5.5, height = 4)
grade3.results %>% plot_contrast_results(grade1.results, cutoff.from = 0.5, cutoff.to = 0.5, trim = 1)
dev.off()


# Reduction in important interactions

imp.inter <- rbind(grade1.results$importances.aggregated %>% 
                     filter(Importance >= 0.5) %>% 
                     group_by(view) %>% 
                     summarize(value = n()) %>% 
                     mutate(grade = "grade 1"),
                   grade3.results$importances.aggregated %>% 
                     filter(Importance >= 0.5) %>% 
                     group_by(view) %>% 
                     summarize(value = n()) %>% 
                     mutate(grade = "grade 3"))

ggplot(imp.inter, aes(x = grade, y = value, group = view, color = view)) + 
  geom_line() + geom_point() + 
  scale_color_brewer(palette = "Set2") +
  xlab("Grade") + ylab("Important interactions") +
  ylim(40,100) +
  theme_classic()

# Permutation analysis

perm.results <- collect_results(seq(10) %>%
  map(~ (list.dirs(paste0("results/imc_small_ridge_perm", .x, "/imc_bc_optim_zoi"))[-1])) %>%
  unlist())

generate_perm_plots <- function(global, meta) {
  ymin <- min(global$value)
  ymax <- max(global$value)

  # global comparison
  global.plot <- ggplot(global, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Global") +
    theme_classic()

  # grades

  grade1 <- global %>% filter(sample %in% (meta %>% filter(Grade == 1) %>% pull(`Sample ID`)))
  grade1.plot <- ggplot(grade1, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 1") +
    theme_classic()

  grade2 <- global %>% filter(sample %in% (meta %>% filter(Grade == 2) %>% pull(`Sample ID`)))
  grade2.plot <- ggplot(grade2, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 2") +
    theme_classic()

  grade3 <- global %>% filter(sample %in% (meta %>% filter(Grade == 3) %>% pull(`Sample ID`)))
  grade3.plot <- ggplot(grade3, aes(x = set, y = value)) +
    geom_violin(aes(fill = set), scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette = "Set2") +
    ylim(ymin, ymax) +
    ggtitle("Grade 3") +
    theme_classic()


  # tests
  
  list(global, grade1, grade2, grade3) %>% map_dbl(function(data) {
    values <- data %>%
      group_by(set) %>%
      group_split() %>%
      map((. %>% pull(value)))
    wilcox.test(values[[1]], values[[2]], alternative = "greater")$p.value %>% log10()
  }) %>% print

  plot_grid(global.plot, grade1.plot, grade2.plot, grade3.plot)
}


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


# difference in gain

global.perf <- bind_rows(
  tibble(set = "original", bc.results$improvements %>% filter(measure == "gain.R2") %>%
    select(sample, value)),
  tibble(set = "permutation", perm.results$improvements %>% filter(measure == "gain.R2") %>%
    select(sample, value))
) %>% mutate(sample = str_extract(sample, "[ABC][a-zA-Z0-9]+"))
generate_perm_plots(global.perf, meta)

# plot coefficient fractions

global.intra  <- get_global_fractions(bc.results, perm.results, "intra")
generate_perm_plots(global.intra, meta)

global.juxta <- get_global_fractions(bc.results, perm.results, "juxta")
generate_perm_plots(global.juxta, meta)

global.para <- get_global_fractions(bc.results, perm.results, "para")
generate_perm_plots(global.para, meta)


frac.significant <- bind_rows(
  tibble(
    set = "original",
    bc.results$contributions %>%
      filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
      group_by(view) %>% summarize(fraction = sum(value < 0.05) / n(), .groups = "drop")
  ),
  tibble(
    set = "permutation",
    perm.results$contributions %>%
      filter(str_starts(view, "p\\."), !str_detect(view, "intercept")) %>%
      group_by(view) %>% summarize(fraction = sum(value < 0.05) / n(), .groups = "drop")
  )
)

ggplot(frac.significant, aes(x = view, y = fraction)) +
  geom_point(aes(color = set), size = 2) +
  ylim(0, 1) +
  theme_classic()



# Singnature analysis

signature <- extract_signature(bc.results, trim = 1)

signature.pca <- prcomp(signature %>% select(-sample) %>% mutate_all(replace_na, 0))

meta.pca <- left_join(as_tibble(signature.pca$x) %>%
                        mutate(sample = signature %>% pull(sample) %>% str_extract("[ABC]y[0-9]+x[0-9]+")),
                      meta %>% filter(`Sample ID` %in% (signature %>% pull(sample) %>% str_extract("[ABC]y[0-9]+x[0-9]+"))),
            by = c("sample" = "Sample ID"))  %>% 
  mutate(Grade = as.factor(Grade))

ggplot(meta.pca %>% filter(!is.na(Grade)), aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Grade), size = 3) + scale_color_brewer(palette = "Set2") + theme_classic()

pdf(file = "plots/imc_small/small_pca_r2_var.pdf", width = 6, height = 5)
fviz_pca_var(signature.pca, col.var = "cos2", repel = TRUE, select.var = list(cos2 = 15),
             gradient.cols = c("#666666","#377EB8", "#E41A1C")) + theme_classic()
dev.off()


# generate importance signature

imp.signature <- extract_signature(bc.results, type = "importance", trim = 1)

impsig.pca.raw <- prcomp(imp.signature %>% select(-sample))

impsig.pca <- data.frame(sample = bc.results$improvements %>% pull(sample) %>% unique %>% 
                      str_extract("[ABC][a-zA-Z0-9]+"), impsig.pca.raw$x)

impmeta.pca <- left_join(meta %>% filter(`Sample ID` %in% (impsig.pca %>% pull(sample))), 
                      impsig.pca, by = c("Sample ID" = "sample")
) %>% mutate(Grade = as.factor(Grade))

ggplot(impmeta.pca %>% filter(!is.na(Grade)), aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Grade), size = 3) + scale_color_brewer(palette = "Set2") + theme_classic()

impsig.pca.raw$rotation

pdf(file = "plots/imc_small/small_pca_imp_var.pdf", width = 9.6, height = 8)
fviz_pca_var(impsig.pca.raw, col.var = "cos2", select.var = list(cos2 = 15),
             gradient.cols = c("#666666","#377EB8", "#E41A1C"), repel = TRUE) + theme_classic()
dev.off()


#raw expression

data <- list.dirs("data/imc_small_breastcancer/", recursive = FALSE)

mean.expr <- data %>% map_dfr(function(d) {
  expr <- read_delim(paste0(d, "/expressions.txt"), delim = " ", col_types = cols())
  colnames(expr) <- make.names(colnames(expr))
  c(image = str_extract(d,"[ABC][a-zA-Z0-9]+"), colMeans(expr))
})

mean.expr.pca <- mean.expr %>% select(-image) %>% mutate_all(as.numeric) %>% prcomp(scale. = FALSE)


mean.expr.join <- left_join(cbind(mean.expr, mean.expr.pca$x), meta, by = c("image" = "Sample ID")) %>% 
  mutate(Grade = as.factor(Grade)) %>% filter(!is.na(Grade))  

ggplot(mean.expr.join, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Grade), size = 3) + scale_color_brewer(palette = "Set2") + theme_classic()

fviz_pca_var(mean.expr.pca, col.var = "cos2", repel = TRUE, col.circle = NA,
             gradient.cols = c("#666666","#377EB8", "#E41A1C")) + theme_classic()
 