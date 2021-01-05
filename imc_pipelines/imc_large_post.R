library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(MISTy)
library(cowplot)
library(ggplot2)
library(uwot)


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

# All grades
ids <- bind_rows(
  clinical.features.bl %>% filter(grade %in% seq(3)) %>%
    select(object.id, grade),
  clinical.features.zh %>% filter(location %in% c("CENTER", "PERIPHERY"), grade %in% seq(3)) %>%
    mutate(grade = as.numeric(grade)) %>%
    select(object.id, grade)
)

allg.results <- collect_results(paste0("results/imc_large_optim/", ids %>% pull(object.id)))

# Grade PCA

r2.signature <- allg.results$improvements %>% filter(str_ends(measure,"R2")) %>% 
  unite("Feature", c(target, measure)) %>% group_by(image) %>%
  pivot_wider(names_from = "Feature", values_from = "value") %>%
  ungroup()

r2.signature.pca <- prcomp(r2.signature %>% select(-image) %>% mutate_all(replace_na, 0), scale. = TRUE)

ggplot(cbind(ids, r2.signature.pca$x), aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = as.factor(grade)), size = 1)  + labs(color = "Grade") + theme_classic()
ggsave("plots/imc_large/pca_grades.pdf")

# Plots by grade
ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 1) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 1) %>%
    pull(object.id)
)

grade1.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/grade_1.pdf", width = 6.5, height = 5)
plot_bundle(grade1.results)
dev.off()

ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 2) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 2) %>%
    pull(object.id)
)

grade2.results <- collect_results(paste0("results/imc_large_optim/", ids))
pdf("plots/imc_large/grade_2.pdf", width = 6.5, height = 5)
plot_bundle(grade2.results)
dev.off()

ids <- c(
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), grade == 3) %>%
    pull(object.id),
  clinical.features.bl %>%
    filter(grade == 3) %>%
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


# Grade 3 importance signature pca/umap
signature <- function(misty.results, view = "para") {
  misty.results$importances %>% map_dfr(~ .x[[view]] %>%
    pivot_longer(names_to = "Target", values_to = "Importance", -Predictor) %>%
    unite("Feature", Predictor, Target) %>%
    filter(!is.na(Importance)) %>%
    pivot_wider(names_from = Feature, values_from = Importance))
}

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

cltype.umap <- umap(cbind(sig_intra, sig_juxta, sig_para), spread = 0.5)
cltype.pca <- prcomp(cbind(sig_intra, sig_juxta, sig_para))
# cltype.umap <- umap(sig_intra*1, spread = 0.5)
colnames(cltype.umap) <- c("U1", "U2")
ggplot(left_join(bind_cols(id = ids %>% pull(object.id), cltype.umap), ids, by = c("id" = "object.id")), aes(x = U1, y = U2, color = clinical_type)) +
  geom_point(size = 2) +
  theme_classic()
ggsave("plots/imc_large/umap_g3_ctype.pdf")
ggplot(left_join(bind_cols(id = ids %>% pull(object.id), cltype.pca$x), ids, by = c("id" = "object.id")), aes(x = PC1, y = PC2, color = clinical_type)) +
  geom_point(size = 2) +
  theme_classic()
ggsave("plots/imc_large/pca_g3_ctype.pdf")

# grade 3 r2 signature pca
r2.signature <- cltype.results$improvements %>%
  filter(str_ends(measure, "R2")) %>%
  unite("Feature", c(target, measure)) %>%
  group_by(image) %>%
  pivot_wider(names_from = "Feature", values_from = "value") %>%
  ungroup()

r2.signature.pca <- prcomp(r2.signature %>% select(-image) %>% mutate_all(replace_na, 0), scale. = TRUE)

ggplot(cbind(ids, r2.signature.pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = clinical_type), size = 1) +
  theme_classic()


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

hrpositive.results <- collect_results(paste0("results/imc_large_optim/", ids))

pdf("plots/imc_large/g3_subtype_HR+.pdf", width = 6.5, height = 5)
plot_bundle(hrpositive.results)
dev.off()

ids <- c(
  clinical.features.bl %>%
    filter(str_detect(clinical_type, "(HR-|TripleNeg)"), grade == 3) %>%
    pull(object.id),
  clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), str_detect(clinical_type, "(HR-|TripleNeg)"), grade == 3) %>%
    pull(object.id)
)

hrnegative.results <- collect_results(paste0("results/imc_large_optim/", ids))

pdf("plots/imc_large/g3_subtype_HR-.pdf", width = 6.5, height = 5)
plot_bundle(hrnegative.results)
dev.off()

pdf("plots/imc_large/contrast_HR+_HR-.pdf", width = 6.5, height = 5)
plot_contrast_results(hrpositive.results, hrnegative.results, cutoff.from = 0.5, cutoff.to = 0.5)
dev.off()

pdf("plots/imc_large/contrast_HR-_HR+.pdf", width = 6.5, height = 5)
plot_contrast_results(hrnegative.results, hrpositive.results, cutoff.from = 0.5, cutoff.to = 0.5)
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
    summarise(Importance = mean(Importance), .groups = "drop_last")
  # summarise(Importance = max(0, Importance), .groups = "drop_last")

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
surv.cor <- c("HR+HER2-", "TripleNeg", "HR+HER2+", "HR-HER2+") %>% map(function(ctype) {
  message(ctype)
  ids.table <- clinical.features.bl %>%
    filter(clinical_type == ctype, grade == 3, OSmonth > 0, Patientstatus == "alive") %>%
    select(object.id, OSmonth) %>%
    rename(feature = OSmonth)

  ids <- ids.table %>% pull(object.id)
  surv.results <- collect_results(paste0("results/imc_large_optim/", ids))

  cor.analysis(ids.table, surv.results, "para")
})

# TODO: plots for correlation analysis


# Playground

# misty.results <- collect_results(paste0("results/imc_large_optim/", c(clinical.features.bl %>%
#                                           pull(object.id), clinical.features.zh %>% filter(location %in% c("CENTER", "PERIPHERY")) %>% pull(object.id))))
#
# # complete set
#
# pdf("complete.pdf", width = 6.5, height = 5)
#
# plot_bundle(misty.results)
#
# dev.off()


# seq(3) %>% walk(function(g){
#   ids <- c(clinical.features.zh %>%
#              #filter(location %in% c("CENTER"), grade == g) %>%
#              ## all locations for consistency
#              filter(location %in% c("CENTER", "PERIPHERY"), grade == g) %>%
#              pull(object.id),
#            clinical.features.bl %>%
#              filter(grade == g) %>%
#              pull(object.id)
#   )
#
#   grade.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
#   pdf(paste0("plots/imc_large/grade_",g,".pdf"), width = 6.5, height = 5)
#   plot_bundle(grade.results)
#   dev.off()
# })

# Plots by location
## Deprecated

# c("CENTER", "PERIPHERY", "[]") %>% walk(function(l){
#   ids <- clinical.features.zh %>%
#     filter(location == l) %>%
#     pull(object.id)
#
#   location.results <- collect_results(paste0("../results/imc_large_optim/", ids))
#
#   pdf(paste0("location_", l, ".pdf"), width = 6.5, height = 5)
#   plot_bundle(location.results)
#   dev.off()
#
# })
#
# ids <- clinical.features.zh %>%
#   filter(location == "[]") %>%
#   pull(object.id)
# normal.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
# ids <- clinical.features.zh %>%
#   filter(location == "PERIPHERY") %>%
#   pull(object.id)
# periphery.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
# ids <- clinical.features.zh %>%
#   filter(location == "CENTER") %>%
#   pull(object.id)
# center.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
# plot_contrast_results(center.results, periphery.results, cutoff.from = 0.5, cutoff.to = 0.5)
# plot_contrast_results(periphery.results, normal.results, cutoff.from = 0.5, cutoff.to = 0.5)
#
#
# # Plots by TNM_N
#
# seq(0,3) %>% walk(function(ndeg){
#   ids <- c(clinical.features.zh %>%
#     filter(location %in% c("CENTER", "PERIPHERY"), str_detect(PTNM_N, as.character(ndeg))) %>%
#     pull(object.id),
#     clinical.features.bl %>%
#       filter(str_detect(PTNM_N, as.character(ndeg))) %>%
#       pull(object.id)
#   )
#
#   ndeg.results <- collect_results(paste0("../results/imc_large_optim/", ids))
#
#   pdf(paste0("ndeg_", ndeg, ".pdf"), width = 6.5, height = 5)
#   plot_bundle(ndeg.results)
#   dev.off()
#
# })


# Interexperiment contrasts
# Grade
# no - plus baseline expressions, structure loss vs interaction

# expression stats
# base_expression_stats <- function(ids){
#   paths <- paste0("../Spatial/spatial_data/JacksonFischer_Collaborators/single_images/ ", ids, ".csv")
#
#   panel <- read_csv("data/imc_large_breastcancer/Basel_Zuri_StainingPanel.csv") %>%
#     select(Target, FullStack) %>% filter(complete.cases(.)) %>% distinct()
#
#   paths %>% map_dfr(function(path){
#     pid <- str_extract(path, "[0-9]+")
#     message(paste0("Processing ", pid))
#
#     data <- read_csv(path)
#     if(nrow(data) < 1000) return()
#
#     expr <- data %>% select(contains("Intensity_MeanIntensity_FullStack"))
#
#     # translate channel to marker
#     markers <- tibble(channel = colnames(expr) %>% str_extract("[0-9]+") %>% as.numeric) %>%
#       left_join(panel %>% slice(-c(41, 46, 47)), by = c("channel" = "FullStack"))
#
#     # cleanup both markers and expr
#     to.remove <- which(markers$channel < 9 | is.na(markers$Target) |
#                          markers$channel > 47 | markers$channel %in% c(12, 25, 28, 38))
#     expr <- expr %>% select(-all_of(to.remove))
#     colnames(expr) <- markers %>% slice(-to.remove) %>% pull(Target) %>% make.names
#
#     # get geometry
#     pos <- data %>% select(Location_Center_X, Location_Center_Y)
#
#     # filter nans
#     pos.complete <- which(complete.cases(pos))
#
#     expr <- expr %>% slice(pos.complete)
#
#     expr
#   })
# }
# bes <- base_expression_stats(ids) %>% pivot_longer(cols = everything(), names_to = "Marker", values_to = "Expression")
# ggplot(bes, aes(x = Marker, y = Expression)) + geom_boxplot(outlier.shape = NA) + ylim(c(0,75)) + theme_classic()



# anova.analysis <- function(ids.table, results,  view = "para", plottop = 4){
#   agg.pair.feat <- ids.table %>% pull(object.id) %>% imap_dfr(function(id, ind){
#     feat <- ids.table %>% pluck("feature", ind)
#     results$importances[[paste0("results/imc_large_optim/", id)]][[view]] %>%
#       pivot_longer(names_to = "Target", values_to = "Importance", -Predictor) %>%
#       filter(!is.na(Importance)) %>% mutate(Feature = feat, Survivability = ids.table %>% pluck("OSmonth", ind))
#   })
#
#   toplot <- agg.pair.feat
#
#   toreturn <- agg.pair.feat %>% mutate(Feature = as.factor(Feature)) %>% group_by(Predictor, Target) %>%
#     #do(anova = TukeyHSD(aov(Importance~Feature, .))$Feature %>% data.frame %>% tibble::rownames_to_column("comb") %>% select(comb, p.adj)) %>%
#     do(anova = unlist(summary(aov(Importance ~ Feature, data = .)))["Pr(>F)1"]) %>%
#     unnest(cols = c(anova)) %>% mutate(q.anova = p.adjust(anova, method = "BH"))
# }



# predictor-target pairs vs overall survivability

# surv.cor <- c("HR+HER2-", "TripleNeg", "HR+HER2+", "HR-HER2+") %>% map(function(ctype){
#   message(ctype)
#   ids.table <- clinical.features.bl %>% filter(clinical_type == ctype, OSmonth > 0) %>% select(object.id, OSmonth) %>% rename(feature = OSmonth)
#
#   ids <- ids.table %>% pull(object.id)
#   surv.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
#   cor.analysis(ids.table, surv.results, "para")
# })
#
# ids.table <- clinical.features.bl %>%
#   filter(OSmonth > 0, !is.na(clinical_type), str_detect(Patientstatus, "alive")) %>%
#   select(object.id, OSmonth, clinical_type)
#
# ids <- ids.table %>% pull(object.id)
# surv.results <- collect_results(paste0("results/imc_large_optim/", ids))
#
#
# cor.results <- cor.analysis(ids.table %>% rename(feature = OSmonth), surv.results, "para")
#
# anova.results <- anova.analysis(ids.table %>% filter(!is.na(feature)), surv.results, "para")
#
#
# featbind <- clinical.features.bl %>%
#   filter(!is.na(clinical_type), str_detect(Patientstatus, "alive")) %>% select(object.id, grade, clinical_type, OSmonth, Patientstatus)
#
# surv.results <- collect_results(paste0("results/imc_large_optim/", featbind %>% pull(object.id)))
# sig <- signature(surv.results)
# sig.filtered <- sig %>% mutate_all(~max(0, .x))
#
# sig.umap <- uwot::umap(sig.filtered, init="normlaplacian", spread = 0.5)
# colnames(sig.umap) <- c("U1", "U2")
# ggplot(left_join(cbind(id = featbind %>% pull(object.id), sig.umap) %>% as_tibble, featbind, by =c("id" = "object.id")) %>% filter(clinical_type == "HR+HER2-"), aes(x = U1, y=U2, color = OSmonth)) + geom_point(size = 3) + theme_classic()
#
#
# signature.gain <- surv.results$improvements %>% filter(str_ends(measure,"R2")) %>%
#   unite("Feature", c(target, measure)) %>% group_by(image) %>%
#   pivot_wider(names_from = "Feature", values_from = "value") %>%
#   mutate(image = str_extract(image, "[0-9]+")) %>% ungroup()
#
# #add also zurich cohort data?
# signature.umap <- uwot::umap(signature.gain %>% select(-image) %>% mutate_all(replace_na, 0), spread = 0.3, min_dist = 1e-5)
# colnames(signature.umap) <- c("U1", "U2")
# ggplot(left_join(cbind(id = featbind %>% pull(object.id), signature.umap) %>% as_tibble, featbind, by =c("id" = "object.id")), aes(x = U1, y=U2)) + geom_point(size = 2) + theme_classic()


# perform correlation stability validation by stratified sampling
# samp.res <- seq(10) %>% map(~
#   cor.analysis(ids.table %>% filter(!is.na(clinical_type)) %>%
#                  group_by(clinical_type) %>%
#                  slice_sample(n=20) %>% ungroup %>%
#                  rename(feature = OSmonth), surv.results, "para") %>%
#     arrange(Predictor, Target)
# )
#
# samp.res[[1]] %>% select(Predictor, Target) %>%
#   tibble::add_column(p.cor = samp.res %>% map_dfc(~.x$p.cor) %>% apply(1, function(x) length(x)/sum(1/x))) %>%
#   mutate(q.cor = p.adjust(p.cor)) %>% arrange(q.cor)




# MISTy to MinHash signature
# generate column index permutation
# reorder the columns based on the index permutation
# find the first true value which.max(sig[n, ])

# minhash_misty <- function(misty.signature, mhsig.length = 50, seed = 1){
#   set.seed(seed)
#
#   mhsig <- seq_len(mhsig.length) %>% map_dfc(function(n){
#     corder <- sample.int(ncol(misty.signature))
#     apply(misty.signature[,corder], 1, which.max)
#   }) %>% `colnames<-`(paste0("MH", seq_len(mhsig.length)))
#
#   mhsig
# }
# cltype.umap <- uwot::umap(cbind(minhash_misty(sig_juxta, mhsig.length = 100), minhash_misty(sig_para, mhsig.length = 100)))
# cltype.umap <- uwot::umap(minhash_misty(cbind(sig_juxta, sig_para), mhsig.length = 500, seed=42))
