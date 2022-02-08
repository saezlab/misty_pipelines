library(tidyverse)
library(mistyR)
library(cowplot)
library(factoextra)
library(gridGraphics)

future::plan(future::multisession)

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

# All grades
ids <- bind_rows(
  clinical.features.bl %>% filter(grade %in% seq(3), !is.na(clinical_type)) %>%
    select(object.id, grade, clinical_type),
  clinical.features.zh %>% filter(location %in% c("CENTER", "PERIPHERY"), grade %in% seq(3), !is.na(clinical_type)) %>%
    mutate(grade = as.numeric(grade)) %>%
    select(object.id, grade, clinical_type)
)

allg.results <- collect_results(paste0("results/imc_large_optim/", ids %>% pull(object.id)))


imp.signature <- extract_signature(allg.results, type = "importance", trim = 1)

imp.signature.pca <- prcomp(imp.signature %>% select(-sample))

f5c <- ggplot(
  left_join(bind_cols(id = ids %>% pull(object.id), imp.signature.pca$x), ids, by = c("id" = "object.id")),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(grade)), size = 1) +
  labs(color = "Grade") +
  theme_classic()


f5d <- ggplot(
  left_join(bind_cols(id = ids %>% pull(object.id), imp.signature.pca$x), ids, by = c("id" = "object.id")),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = clinical_type), size = 1) +
  labs(color = "Clinical type") +
  theme_classic()


# For supplementary
allg.results %>%
  plot_improvement_stats() %>%
  plot_improvement_stats("intra.R2") %>%
  plot_view_contributions(trim = 1)

fviz_pca_var(imp.signature.pca,
  col.var = "cos2", select.var = list(cos2 = 15), repel = TRUE,
  gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
) + theme_classic()


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

pdf("plots/imc_large/community1.pdf")
grade1.results %>% plot_interaction_communities(view = "juxta", cutoff = 0.5)
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

pdf("plots/imc_large/community2.pdf")
grade2.results %>% plot_interaction_communities(view = "juxta", cutoff = 0.5)
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

pdf("plots/imc_large/community3.pdf")
grade3.results %>% plot_interaction_communities(view = "juxta", cutoff = 0.5)
dev.off()


imp.inter <- rbind(
  grade1.results$importances.aggregated %>%
    filter(Importance >= 0.5) %>%
    group_by(view) %>%
    summarize(value = n()) %>%
    mutate(grade = "grade 1"),
  grade2.results$importances.aggregated %>%
    filter(Importance >= 0.5) %>%
    group_by(view) %>%
    summarize(value = n()) %>%
    mutate(grade = "grade 2"),
  grade3.results$importances.aggregated %>%
    filter(Importance >= 0.5) %>%
    group_by(view) %>%
    summarize(value = n()) %>%
    mutate(grade = "grade 3")
)

f5a <- ggplot(imp.inter, aes(x = grade, y = value, group = view, color = view)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  xlab("Grade") +
  ylab("Important interactions") +
  theme_classic()


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

g3.signature <- extract_signature(grade3.results, type = "importance", trim = 1)

cltype.pca <- prcomp(g3.signature %>% select(-sample), scale. = TRUE)

ggplot(left_join(bind_cols(id = ids %>% pull(object.id), cltype.pca$x), ids, by = c("id" = "object.id")), aes(x = PC1, y = PC2, color = clinical_type)) +
  geom_point(size = 2) +
  theme_classic()

fviz_pca_var(cltype.pca,
  col.var = "cos2", select.var = list(cos2 = 12), repel = TRUE,
  gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
) + theme_classic()


# Correlation analysis of importance to survivability
cor.analysis <- function(ids.table, results, view = "para", plottop = 4, trim = 1) {
  agg.pair.feat <- ids.table %>%
    pull(object.id) %>%
    imap_dfr(function(id, ind) {
      feat <- ids.table %>% pluck("feature", ind)
      targets <- results$improvements.stats %>%
        filter(measure == "gain.R2", mean >= trim) %>%
        pull(target)
      results$importances %>%
        filter(
          sample == paste0(getwd(), "/results/imc_large_optim/", id),
          Target %in% targets,
          view == !!view,
          !is.na(Importance)
        ) %>%
        mutate(Feature = feat)
    })

  toplot <- agg.pair.feat %>%
    group_by(Predictor, Target, Feature) %>%
    summarise(Importance = mean(Importance), .groups = "drop_last") %>%
    filter(sum(Importance >= 0) >= 0.3 * n()) # %>%

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


# All
ids.table <- clinical.features.bl %>%
  filter(grade == 3, OSmonth > 0, Patientstatus == "alive") %>%
  select(object.id, OSmonth) %>%
  rename(feature = OSmonth)

ids <- ids.table %>% pull(object.id)
surv.results <- collect_results(paste0("results/imc_large_optim/", ids))

all.surv.cor <- bind_rows(
  cor.analysis(ids.table, surv.results, "intra") %>% mutate(view = "intra"),
  cor.analysis(ids.table, surv.results, "juxta") %>% mutate(view = "juxta"),
  cor.analysis(ids.table, surv.results, "para") %>% mutate(view = "para")
)

# for grade 3 tumors
subtypes <- c("HR+HER2-", "TripleNeg", "HR+HER2+", "HR-HER2+")

surv.cor <- subtypes %>%
  map(function(ctype) {
    message(ctype)
    ids.table <- clinical.features.bl %>%
      filter(clinical_type == ctype, grade == 3, OSmonth > 0, Patientstatus == "alive") %>%
      select(object.id, OSmonth) %>%
      rename(feature = OSmonth)

    ids <- ids.table %>% pull(object.id)
    surv.results <- collect_results(paste0("results/imc_large_optim/", ids))

    bind_rows(
      cor.analysis(ids.table, surv.results, "intra") %>% mutate(view = "intra"),
      cor.analysis(ids.table, surv.results, "juxta") %>% mutate(view = "juxta"),
      cor.analysis(ids.table, surv.results, "para") %>% mutate(view = "para")
    )
  }) %>%
  `names<-`(subtypes)

all.surv.cor %>%
  filter(p.cor <= 0.05) %>%
  arrange(desc(abs(cor))) %>%
  select(-q.cor) %>%
  write_csv("top_cor.txt")

surv.cor %>%
  map_dfr(~ .x %>%
    filter(p.cor <= 0.05) %>%
    arrange(desc(abs(cor))) %>%
    select(-q.cor), .id = "Subtype") %>%
  write_csv("top_cor_subtypes.txt")



# Kaplan - Meier

km_curve <- function(results, surv, view, predictor, target, cutoff = 1) {
  # collect p-t importance from results
  ptimp <- results$importances %>% 
    filter(view == view, Predictor == predictor, Target == target) %>% 
    mutate(sample = str_extract(sample, "[0-9]+$"))

  group.sum <- left_join(surv,
    tibble(
      object.id = ptimp  %>% pull("sample") %>% as.double(),
      ptimp = ptimp %>% pull("Importance")
    ),
    by = "object.id"
  ) %>%
    mutate(ptimp = ptimp < cutoff) %>%
    group_by(ptimp, feature) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    arrange(feature) %>%
    mutate(pct = 100 - cumsum(n * 100 / sum(n))) %>%
    ungroup()

  dummy <- tribble(
    ~ptimp, ~feature, ~n, ~pct,
    TRUE, 0, 0, 100,
    FALSE, 0, 0, 100
  )

  J <- group.sum %>%
    pull(feature) %>%
    unique()

  hazard <- group.sum %>%
    select(-pct) %>%
    add_row(expand.grid(ptimp = c(TRUE, FALSE), feature = c(0, J), n = 0)) %>%
    group_by(ptimp, feature) %>%
    summarize(n = sum(n), .groups = "drop_last") %>%
    mutate(cn = cumsum(n), N = sum(n) - cn) %>%
    group_split() %>%
    map(~ .x %>% select(-c(ptimp, cn)))

  distribution <- full_join(hazard[[1]], hazard[[2]], by = "feature", suffix = c("i", "j")) %>%
    mutate(
      O = ni + nj, N = Ni + Nj,
      Ei = Ni * O / N, Vi = Ei * ((N - O) / N) * ((N - Ni) / (N - 1)) # ,
      # Ej = Nj*O/N,
      # Vj = Ej*((N-O)/N)*((N-Nj)/(N-1))
    ) %>%
    drop_na()

  Zi <- sum(distribution$ni - distribution$Ei) / sqrt(sum(distribution$Vi))

  pval <- format(round(pnorm(-abs(Zi)), 3), nsmall = 3)

  ggplot(group.sum %>% add_row(dummy), aes(x = feature, y = pct, color = ptimp)) +
    geom_step() +
    geom_point(shape = 3) +
    labs(
      x = "Time (Months)", y = "Percent survival", color = "P-T interaction",
      title = paste(view, predictor, target, pval)
    ) +
    theme_classic()
}

top3 <- all.surv.cor %>%
  filter(p.cor <= 0.05) %>%
  arrange(desc(abs(cor))) %>%
  group_by(view) %>% slice(1)


kms <- top3 %>% pmap(function(...){
  current <- tibble(...)
  km_curve(surv.results, ids.table, current$view, current$Predictor, current$Target)
})

f5e <- kms[[1]]
f5f <- kms[[2]]
f5g <- kms[[3]]

# Figure 5
pdf("plots/imc_large/Figure5.pdf", width = 13.2, height = 16.6)
plot_grid(f5a, grid::nullGrob(), grid::nullGrob(), grid::nullGrob(), 
          plot_grid(f5c, f5d, ncol = 1, labels = c("C", "D")), f5e, f5f, f5g, 
          ncol = 2, byrow = FALSE, rel_heights = c(2,1,1,1),
          labels = c("A", "", "B", "E", "", "F", "", "G"))
dev.off()


# Supplementary
 

