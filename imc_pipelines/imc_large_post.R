library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(MISTy)


wi <- read_csv("../data/imc_large_breastcancer/Basel_Zuri_WholeImage.csv")
bm <- read_csv("../data/imc_large_breastcancer/Basel_PatientMetadata.csv")
zm <- read_csv("../data/imc_large_breastcancer/Zuri_PatientMetadata.csv")

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

plot_bundle <- function(misty.results){
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

misty.results <- collect_results(paste0("../results/imc_large_optim/", c(clinical.features.bl %>% 
                                          pull(object.id), clinical.features.zh %>% filter(location %in% c("CENTER", "PERIPHERY")) %>% pull(object.id))))

# complete set

pdf("complete.pdf", width = 6.5, height = 5)

plot_bundle(misty.results)

dev.off()


# Plots by grade
seq(3) %>% walk(function(g){
  ids <- c(clinical.features.zh %>%
             filter(location %in% c("CENTER", "PERIPHERY"), grade == g) %>%
             pull(object.id),
           clinical.features.bl %>%
             filter(grade == g) %>%
             pull(object.id)
  )
  
  grade.results <- collect_results(paste0("../results/imc_large_optim/", ids))
  
  pdf(paste0("grade_",g,".pdf"), width = 6.5, height = 5)
  plot_bundle(grade.results)
  dev.off()
  
})

# Plots by location

c("CENTER", "PERIPHERY", "[]") %>% walk(function(l){
  ids <- clinical.features.zh %>%
    filter(location == l) %>%
    pull(object.id)
  
  location.results <- collect_results(paste0("../results/imc_large_optim/", ids))
  
  pdf(paste0("location_", l, ".pdf"), width = 6.5, height = 5)
  plot_bundle(location.results)
  dev.off()
  
})



# Plots by TNM_N

seq(0,3) %>% walk(function(ndeg){
  ids <- c(clinical.features.zh %>%
    filter(location %in% c("CENTER", "PERIPHERY"), str_detect(PTNM_N, as.character(ndeg))) %>%
    pull(object.id),
    clinical.features.bl %>%
      filter(str_detect(PTNM_N, as.character(ndeg))) %>%
      pull(object.id)
  )
  
  ndeg.results <- collect_results(paste0("../results/imc_large_optim/", ids))
  
  pdf(paste0("ndeg_", ndeg, ".pdf"), width = 6.5, height = 5)
  plot_bundle(ndeg.results)
  dev.off()

})

# Plots by survivability

surv1 <- c(0,3,10,22,39)
surv2 <- c(3,10,22,39,158)

surv1 %>% walk2(surv2, function(from, to){
  ids <- clinical.features.bl %>% mutate(surv = OSmonth - DFSmonth) %>%
    filter(surv > from, surv <=to) %>%
    pull(object.id)
  
  surv.results <- collect_results(paste0("../results/imc_large_optim/", ids))
  
  pdf(paste0("survivability_", from, "_", to, ".pdf"), width = 6.5, height = 5)
  plot_bundle(surv.results)
  dev.off()
})
