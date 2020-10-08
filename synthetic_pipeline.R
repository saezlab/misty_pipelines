library(future)
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(MISTy)
library(scales)
library(cowplot)
library(PRROC)


# Run MISTy
data <- list.dirs("data/synthetic", recursive = FALSE)
plan(multiprocess, workers = 4)

l <- 10


data %>% walk(function(d) {
  all <- read_csv(paste0(d, "/random1_position_expression_real_cells.csv")) %>%
    select(-type)
  expr <- all %>% select(-row, -col)
  pos <- all %>% select(row, col)

  views <- create_initial_view(expr) %>% add_paraview(pos, l^2)

  run_misty(views, results.folder = paste0(
    "results/synthetic/",
    str_extract(d, "synthetic[0-9]+"), "/"
  ))
})


# Analysis
results <- list.dirs("results/synthetic", recursive = FALSE)
misty.results <- collect_results(results)


# Plot contributions and improvement
misty.results %>% plot_view_contributions() %>%
  plot_improvement_stats()

# Plot aggregated importances
misty.results %>% plot_interaction_heatmap("intra", 0.5) %>% 
  plot_interaction_heatmap("para.100", 0.5) %>% 
  plot_contrast_heatmap("intra", "para.100")

misty.results %>% plot_interaction_communities(view = "intra", cutoff = 0.5) %>%
  plot_interaction_communities(view = "para.100", cutoff = 0.5)

true.connections <- read_csv("data/synthetic/true_connections.csv")

# Plot true connections
ggtrue.intra <- ggplot(true.connections %>% filter(view == "intra")) + 
  geom_tile(aes(x = node1, y = node2, fill = as_factor(present))) +
  scale_fill_discrete(type = c("white", muted("blue"))) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("true.intra")

ggtrue.para <- ggplot(true.connections %>% filter(view == "para")) + 
  geom_tile(aes(x = node1, y = node2, fill = as_factor(present))) +
  scale_fill_discrete(type = c("white", muted("blue"))) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("true.para")

plot_grid(ggtrue.intra, ggtrue.para)



# Distributions of AUROC and AUPRC

intra.roc <- seq(10) %>% map(function(i){
  tidy.intra <- misty.results$importances[[i]][["intra"]] %>% 
    pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
  joined <- true.connections %>% filter(view == "intra") %>% 
    left_join(tidy.intra, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
    filter(!is.na(present), !is.na(Prediction))
  roc.curve(joined %>% pull(Prediction), weights.class0 = joined %>% pull(present))$auc
})

intra.pr <- seq(10) %>% map(function(i){
  tidy.intra <- misty.results$importances[[i]][["intra"]] %>% 
    pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
  joined <- true.connections %>% filter(view == "intra") %>% 
    left_join(tidy.intra, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
    filter(!is.na(present), !is.na(Prediction))
  pr.curve(joined %>% pull(Prediction), weights.class0 = joined %>% pull(present))$auc.integral
})

para.roc <- seq(10) %>% map(function(i){
  tidy.para <- misty.results$importances[[i]][["para.100"]] %>% 
    pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
  joined <- true.connections %>% filter(view == "para") %>% 
    left_join(tidy.para, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
    filter(!is.na(present), !is.na(Prediction))
  roc.curve(joined %>% pull(Prediction), weights.class0 = joined %>% pull(present))$auc
})

para.pr <- seq(10) %>% map(function(i){
  tidy.para <- misty.results$importances[[i]][["para.100"]] %>% 
    pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
  joined <- true.connections %>% filter(view == "para") %>% 
    left_join(tidy.para, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
    filter(!is.na(present), !is.na(Prediction))
  pr.curve(joined %>% pull(Prediction), weights.class0 = joined %>% pull(present))$auc.integral
})

viodata.roc <- data.frame(intra.roc = intra.roc %>% unlist, para.roc = para.roc %>% unlist) %>% 
  pivot_longer(names_to = "Type", values_to = "AUC")

ggplot(viodata.roc, aes(x=Type, y=AUC)) + geom_violin(aes(fill=Type)) + 
  geom_jitter() + 
  theme_classic() + ylim(0.5,1)

intra.intercept <- sum(true.connections %>% filter(view == "intra", !is.na(present)) %>% pull(present))/
  (true.connections %>% filter(view == "intra", !is.na(present)) %>% nrow())
para.intercept <- sum(true.connections %>% filter(view == "para", !is.na(present)) %>% pull(present))/
  (true.connections %>% filter(view == "para", !is.na(present)) %>% nrow())

viodata.pr <- data.frame(intra.pr = intra.pr %>% unlist, para.pr = para.pr %>% unlist) %>%
  pivot_longer(names_to = "Type", values_to = "AUC")

ggplot(viodata.pr, aes(x=Type, y=AUC)) + geom_violin(aes(fill=Type)) + 
  geom_jitter() + 
  geom_hline(yintercept=intra.intercept, color="#F8766D", linetype="dashed") + 
  geom_hline(yintercept=para.intercept, color="#00BFC4", linetype="dashed") +
  theme_classic() + ylim(0,1)


# Plot aggregated ROC and PR curves

# for iso-F1 and iso-J lines
r <- seq(0,1.2,by=1e-3)
ks <- c(0.1,0.2,0.5,0.8)

tidy.intra.agg <- misty.results$importances.aggregated[["intra"]] %>% 
  pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
joined.intra.agg <- true.connections %>% filter(view == "intra") %>% 
  left_join(tidy.intra.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
  filter(!is.na(present), !is.na(Prediction))

plot(roc.curve(joined.intra.agg %>% pull(Prediction), weights.class0 = joined.intra.agg %>% pull(present), 
               curve = TRUE, rand.compute = TRUE), color="black", rand.plot = TRUE)

ks %>% walk(function(k){
  s <- r + k
  s[s>1] <- NA
  lines(r, s, ylim=c(0,1),xlim=c(0,1), col="gray80")
})

plot(pr.curve(joined.intra.agg %>% pull(Prediction), weights.class0 = joined.intra.agg %>% pull(present), 
              curve = TRUE, rand.compute = TRUE), color="black", rand.plot = TRUE)

ks %>% walk(function(k){
  p <- k*r/(2*r - k)
  p[p<=0 | p > 1.5] <- NA
  lines(r, p, ylim=c(0,1),xlim=c(0,1), col="gray80")
})


tidy.para.agg <- misty.results$importances.aggregated[["para.100"]] %>% 
  pivot_longer(names_to = "Target", values_to = "Prediction", -Predictor)
joined.para.agg <- true.connections %>% filter(view == "para") %>% 
  left_join(tidy.para.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
  filter(!is.na(present), !is.na(Prediction))

plot(roc.curve(joined.para.agg %>% pull(Prediction), weights.class0 = joined.para.agg %>% pull(present), 
               curve = TRUE, rand.compute = TRUE), color="black", rand.plot = TRUE)

ks %>% walk(function(k){
  s <- r + k
  s[s>1] <- NA
  lines(r, s, ylim=c(0,1),xlim=c(0,1), col="gray80")
})

plot(pr.curve(joined.para.agg %>% pull(Prediction), weights.class0 = joined.para.agg %>% pull(present), 
              curve = TRUE, rand.compute = TRUE), color="black", rand.plot = TRUE)

ks %>% walk(function(k){
  p <- k*r/(2*r - k)
  p[p<=0 | p > 1.5] <- NA
  lines(r, p,ylim=c(0,1),xlim=c(0,1), col="gray80")
})