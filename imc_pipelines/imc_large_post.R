library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(MISTy)
library(ggplot2)
library(cowplot)
library(factoextra)


wi <- read_csv("data/imc_large_breastcancer/Basel_Zuri_WholeImage.csv")
bm <- read_csv("data/imc_large_breastcancer/Basel_PatientMetadata.csv")
zm <- read_csv("data/imc_large_breastcancer/Zuri_PatientMetadata.csv")

# Basel mapping
object.id.bl <- seq(nrow(bm)) %>% map_dbl(function(i) {
  file.pattern <- paste0(
    str_extract(
      bm %>% pluck("FileName_FullStack", i),
      "(.+?_){7}"
    ), "[0-9]+_",
    str_extract(bm %>% pluck("core", i), "X.+")
  )

  wi %>%
    pluck("ImageNumber", str_which(wi %>% pull(FileName_FullStack), file.pattern))
})

# Zurich mapping
object.id.zh <- seq(nrow(zm)) %>% map_dbl(function(i) {
  file.pattern <- paste0(
    str_extract(
      zm %>% pluck("FileName_FullStack", i),
      "(.+?_){7}"
    ), "B[0-9]{2}\\.[0-9]+_",
    str_extract(zm %>% pluck("FileName_FullStack", i), "[ABC]y[0-9]+x[0-9]+_[0-9]+_")
  )

  wi %>%
    pluck("ImageNumber", str_which(wi %>% pull(FileName_FullStack), file.pattern))
})



# pesky 307
clinical.features.bl <- bm %>%
  mutate(object.id = object.id.bl, .after = 1) %>%
  filter(Count_Cells >= 1000, diseasestatus == "tumor", object.id != 307)

clinical.features.zh <- zm %>%
  mutate(object.id = object.id.zh, .after = 1) %>%
  filter(Count_Cells >= 1000)

plot_bundle <- function(misty.results) {
  misty.results %>%
    plot_improvement_stats() %>%
    plot_view_contributions() %>%
    plot_interaction_heatmap("intra", 0.5) %>%
    plot_interaction_heatmap("juxta", 0.5) %>%
    plot_interaction_heatmap("para", 0.5) %>%
    plot_interaction_communities("intra", 0.5) %>%
    plot_interaction_communities("juxta", 0.5) %>%
    plot_interaction_communities("para", 0.5) %>%
    plot_contrast_heatmap("intra", "juxta", 0.5) %>%
    plot_contrast_heatmap("intra", "para", 0.5)
}

signature <- function(misty.results, view = "para") {
  misty.results$importances %>% map_dfr(~ .x[[view]] %>%
                                          pivot_longer(names_to = "Target", values_to = "Importance", -Predictor) %>%
                                          unite("Feature", Predictor, Target) %>%
                                          filter(!is.na(Importance)) %>%
                                          pivot_wider(names_from = Feature, values_from = Importance))
}

# All grades
ids <- bind_rows(
  clinical.features.bl %>% filter(grade %in% seq(3), !is.na(clinical_type)) %>%
    select(object.id, grade, clinical_type),
  clinical.features.zh %>% filter(location %in% c("CENTER", "PERIPHERY"), grade %in% seq(3), !is.na(clinical_type)) %>%
    mutate(grade = as.numeric(grade)) %>%
    select(object.id, grade, clinical_type)
)

allg.results <- collect_results(paste0("../results/imc_large_optim/", ids %>% pull(object.id)))
pdf("plots/imc_large/all.pdf", width = 6.5, height = 5)
plot_bundle(allg.results)
dev.off()

# Grade/Clinical type PCA

r2.signature <- allg.results$improvements %>% filter(str_ends(measure,"R2"), !str_ends(measure,"p.R2")) %>% 
  unite("Feature", c(target, measure)) %>% group_by(image) %>%
  pivot_wider(names_from = "Feature", values_from = "value") %>%
  ungroup()

r2.signature.pca <- prcomp(r2.signature %>% select(-image) %>% mutate_all(replace_na, 0), scale. = TRUE)

ggplot(cbind(ids, r2.signature.pca$x), aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = as.factor(grade)), size = 1)  + labs(color = "Grade") + theme_classic()
ggsave("plots/imc_large/pca_r2_grades.pdf", width = 5.5, height = 4)

ggplot(cbind(ids, r2.signature.pca$x), aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = clinical_type), size = 1)  + labs(color = "Clinical type") + theme_classic()
ggsave("plots/imc_large/pca_r2_ctype.pdf", width = 6, height = 4)

fviz_pca_var(r2.signature.pca, col.var = "cos2", select.var = list(cos2 = 10), repel = TRUE, 
             gradient.cols = c("#666666","#377EB8", "#E41A1C"), col.circle = NA) + theme_classic()
ggsave("plots/imc_large/pca_r2_var.pdf", width = 5.5, height = 4)


sig_intra <- signature(allg.results, "intra") # > threshold
sig_juxta <- signature(allg.results, "juxta") # > threshold
sig_para <- signature(allg.results, "para") # > threshold

imp.signature.pca <- prcomp(bind_cols(sig_intra %>% rename_all(~paste0("intra_",.)), 
                               sig_juxta %>% rename_all(~paste0("juxta_",.)), 
                               sig_para %>% rename_all(~paste0("para_",.))))

ggplot(left_join(bind_cols(id = ids %>% pull(object.id), imp.signature.pca$x), ids, by = c("id" = "object.id")), 
       aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(grade)), size = 1) +
  labs(color = "Grade") +
  theme_classic()
ggsave("plots/imc_large/pca_imp_grades.pdf", width = 5.5, height = 4)

ggplot(left_join(bind_cols(id = ids %>% pull(object.id), imp.signature.pca$x), ids, by = c("id" = "object.id")), 
       aes(x = PC1, y = PC2)) +
  geom_point(aes(color = clinical_type), size = 1) +
  labs(color = "Clinical type") +
  theme_classic()
ggsave("plots/imc_large/pca_imp_ctype.pdf", width = 6, height = 4)

fviz_pca_var(imp.signature.pca, col.var = "cos2", select.var = list(cos2 = 10), repel = TRUE, 
             gradient.cols = c("#666666","#377EB8", "#E41A1C"), col.circle = NA) + theme_classic()
ggsave("plots/imc_large/pca_imp_var.pdf", width = 5.5, height = 4)


# Plots by grade
ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 1, !is.na(clinical_type)) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 1, !is.na(clinical_type)) %>%
    pull(object.id)
)

grade1.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/grade_1.pdf", width = 6.5, height = 5)
plot_bundle(grade1.results)
dev.off()

ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 2, !is.na(clinical_type)) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 2, !is.na(clinical_type)) %>%
    pull(object.id)
)

grade2.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/grade_2.pdf", width = 6.5, height = 5)
plot_bundle(grade2.results)
dev.off()

ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 3, !is.na(clinical_type)) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 3, !is.na(clinical_type)) %>%
    pull(object.id)
)

grade3.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/grade_3.pdf", width = 6.5, height = 5)
plot_bundle(grade3.results)
dev.off()

# Grade contrasts
pdf("plots/imc_large/contrast_g3_g2.pdf", width = 6.5, height = 5)
plot_contrast_results(grade3.results, grade2.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

pdf("plots/imc_large/contrast_g2_g1.pdf", width = 6.5, height = 5)
plot_contrast_results(grade2.results, grade1.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

pdf("plots/imc_large/contrast_g3_g1.pdf", width = 6.5, height = 5)
plot_contrast_results(grade3.results, grade1.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

# plots by HR subtype

#HR+
ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), str_detect(clinical_type, "HR+")) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(str_detect(clinical_type, "HR+")) %>%
    pull(object.id)
)

hrp.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/HR+.pdf", width = 6.5, height = 5)
plot_bundle(hrp.results)
dev.off()

ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), !is.na(clinical_type), !str_detect(clinical_type, "HR+")) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(!is.na(clinical_type), !str_detect(clinical_type, "HR+")) %>%
    pull(object.id)
)

hrn.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/HR-.pdf", width = 6.5, height = 5)
plot_bundle(hrn.results)
dev.off()

pdf("plots/imc_large/contrast_HR+_HR-.pdf", width = 6.5, height = 5)
plot_contrast_results(hrp.results, hrn.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

pdf("plots/imc_large/contrast_HR-_HR+.pdf", width = 6.5, height = 5)
plot_contrast_results(hrn.results, hrp.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

# Grade 3 importance signature pca

ids <- bind_rows(
  clinical.features.bl %>%
    filter(!is.na(clinical_type), grade == 3) %>%
    select(object.id, clinical_type),
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), !is.na(clinical_type), grade == 3) %>%
    mutate(grade = as.numeric(grade)) %>%
    select(object.id, clinical_type)
)

cltype.results <- collect_results(paste0("results/imc_large_optim/", ids %>% pull(object.id)))

# threshold <- 0.5
sig_intra <- signature(cltype.results, "intra") # > threshold
sig_juxta <- signature(cltype.results, "juxta") # > threshold
sig_para <- signature(cltype.results, "para") # > threshold

cltype.pca <- prcomp(bind_cols(sig_intra %>% rename_all(~paste0("intra_",.)), 
                               sig_juxta %>% rename_all(~paste0("juxta_",.)), 
                               sig_para %>% rename_all(~paste0("para_",.))))

ggplot(left_join(bind_cols(id = ids %>% pull(object.id), cltype.pca$x), ids, by = c("id" = "object.id")), aes(x = PC1, y = PC2, color = clinical_type)) +
  geom_point(size = 2) +
#  scale_color_brewer(palette = "Set2") +
  theme_classic()
ggsave("plots/imc_large/pca_g3_imp_ctype.pdf", width = 6, height = 4)

fviz_pca_var(cltype.pca, col.var = "cos2", select.var = list(cos2 = 10), repel = TRUE, 
             gradient.cols = c("#666666","#377EB8", "#E41A1C"), col.circle = NA) + theme_classic()
ggsave("plots/imc_large/pca_g3_imp_var.pdf", width = 5.5, height = 4)


# grade 3 r2 signature pca
r2.signature <- cltype.results$improvements %>%
  filter(str_ends(measure, "R2"), !str_ends(measure, "p.R2")) %>%
  unite("Feature", c(target, measure)) %>%
  group_by(image) %>%
  pivot_wider(names_from = "Feature", values_from = "value") %>%
  ungroup()

r2.signature.pca <- prcomp(r2.signature %>% select(-image) %>% mutate_all(replace_na, 0), scale. = TRUE)

ggplot(cbind(ids, r2.signature.pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = clinical_type), size = 2) + 
#  scale_color_brewer(palette = "Set2") +
  theme_classic()

ggsave("plots/imc_large/pca_g3_r2_ctype.pdf", width = 6, height = 4)

fviz_pca_var(r2.signature.pca, col.var = "cos2", select.var = list(cos2 = 10), repel = TRUE, 
             gradient.cols = c("#666666","#377EB8", "#E41A1C"), col.circle = NA) + theme_classic()
ggsave("plots/imc_large/pca_g3_r2_var.pdf", width = 5.5, height = 4)


# Plots by clinical subtype
# c("HR+HER2-", "TripleNeg", "HR+HER2+", "HR-HER2+")

ids <- c(
  clinical.features.bl %>%
    filter(str_detect(clinical_type, "HR+"), grade == 3) %>%
    pull(object.id),
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), str_detect(clinical_type, "HR+"), grade == 3) %>%
    pull(object.id)
)

hrp.results <- collect_results(paste0("results/imc_large_optim/", ids))

pdf("plots/imc_large/g3_subtype_HR+.pdf", width = 6.5, height = 5)
plot_bundle(hrp.results)
dev.off()

ids <- c(
  clinical.features.bl %>%
    filter(str_detect(clinical_type, "(HR-|TripleNeg)"), grade == 3) %>%
    pull(object.id),
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), str_detect(clinical_type, "(HR-|TripleNeg)"), grade == 3) %>%
    pull(object.id)
)

hrn.results <- collect_results(paste0("results/imc_large_optim/", ids))

pdf("plots/imc_large/g3_subtype_HR-.pdf", width = 6.5, height = 5)
plot_bundle(hrn.results)
dev.off()

pdf("plots/imc_large/contrast_g3_HR+_HR-.pdf", width = 6.5, height = 5)
plot_contrast_results(hrp.results, hrn.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

pdf("plots/imc_large/contrast_g3_HR-_HR+.pdf", width = 6.5, height = 5)
plot_contrast_results(hrn.results, hrp.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

# Correlation analysis of importance to survivability
cor.analysis <- function(ids.table, results, view = "para", plottop = 4) {
  agg.pair.feat <- ids.table %>%
    pull(object.id) %>%
    imap_dfr(function(id, ind) {
      feat <- ids.table %>% pluck("feature", ind)
      results$importances[[paste0("results/imc_large_optim/", id)]][[view]] %>%
        pivot_longer(names_to = "Target", values_to = "Importance", -Predictor) %>%
        filter(!is.na(Importance)) %>%
        mutate(Feature = feat)
    })

  toplot <- agg.pair.feat %>%
    group_by(Predictor, Target, Feature) %>%
    summarise(Importance = mean(Importance), .groups = "drop_last") %>% 
    filter(sum(Importance >= 0) >= 0.3*n()) #%>% 
    #filter(n() >= 5)
    #summarise(Importance = max(0, mean(Importance)), .groups = "drop_last")

  toreturn <- toplot %>%
    summarise(
      cor = cor(Importance, Feature, method = "spearman"),
      p.cor = cor.test(Importance, Feature, method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(q.cor = p.adjust(p.cor, method = "BH")) %>%
    arrange(p.cor)

  seq(plottop) %>%
    map(~ toplot %>%
      ungroup() %>%
      filter(
        Predictor == toreturn %>% pluck("Predictor", .x),
        Target == toreturn %>% pluck("Target", .x)
      ) %>%
      ggplot(aes(x = Feature, y = Importance)) +
      geom_point() +
      xlab("Survivability") +
      ggtitle(toreturn[.x, ] %>% paste(collapse = " ")) +
      theme_classic()) %>%
    plot_grid(plotlist = .) %>%
    print()


  return(toreturn)
}

# for grade 3 tumors
subtypes <- c("HR+HER2-", "TripleNeg", "HR+HER2+", "HR-HER2+")

surv.cor <- subtypes %>% map(function(ctype) {
  message(ctype)
  ids.table <- clinical.features.bl %>%
    filter(clinical_type == ctype, grade == 3, OSmonth > 0, Patientstatus == "alive") %>%
    select(object.id, OSmonth) %>%
    rename(feature = OSmonth)

  ids <- ids.table %>% pull(object.id)
  surv.results <- collect_results(paste0("results/imc_large_optim/", ids))

  bind_rows(cor.analysis(ids.table, surv.results, "intra") %>% mutate(view = "intra"),
  cor.analysis(ids.table, surv.results, "juxta") %>% mutate(view = "juxta"),
  cor.analysis(ids.table, surv.results, "para") %>% mutate(view = "para"))
}) %>% `names<-`(subtypes)

# TODO: plots for correlation analysis
