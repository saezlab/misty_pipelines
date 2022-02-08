library(future)
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(mistyR)
library(scales)
library(cowplot)
library(PRROC)
library(ggplot2)
library(tidyr)
source("./large_synthetic_helpers.R")

input_folder = "data/TIGER_tissue_100/"
results_folder_base = "results/TIGER_tissue_100/"
plot_folder_base = "plots/TIGER_tissue_100/"

# Run MISTy
data <- list.dirs(input_folder, recursive = FALSE)
plan(multisession)

l <- 10

view_filter <- function(views, cells) {
  message("Filtering")
  views %>% map(function(view) {
    if (length(view) > 1) {
      view$data <- view$data %>% slice(cells)
    }
    view
  })
}


data %>% walk(function(d) {
  print(paste("processing:",d))
  
  all <- read_csv(paste0(d, "/time100_position_expression_real_cells.csv")) %>%
    select(-starts_with("L"))
  expr <- all %>% select(-row, -col, -type)
  pos <- all %>% select(row, col)
  
  views <- create_initial_view(expr) %>% add_paraview(pos, l)
  
  run_misty(views, results.folder = paste0(results_folder_base,
                                           "all/",
                                           str_extract(d, "IMC[0-9]+"), "/"
  ))
  
  types <- all %>% pull(type)
  
  unique(types) %>% walk(function(t) {
    run_misty(views %>% filter_views(which(types == t)),
              results.folder = paste0(
                results_folder_base, t, "/",
                str_extract(d, "IMC[0-9]+"), "/"
              )
    )
  })
})

# 
# data %>% walk(function(d) {
#   all <- read_csv(paste0(d, "/random1_position_expression_real_cells.csv")) %>%
#     select(-starts_with("lig"))
#   expr <- all %>% select(-row, -col, -type)
#   pos <- all %>% select(row, col)
#   
#   views <- create_initial_view(expr) %>% add_paraview(pos, l)
#   
#   run_misty(views, results.folder = paste0(
#     results_folder,
#     str_extract(d, "_IMC[0-9]+"), "/"
#   ))
# })


ct_results_folder = list(paste0(results_folder_base,"all"),
                         paste0(results_folder_base,"CT1"),
                         paste0(results_folder_base,"CT2"),
                         paste0(results_folder_base,"CT3"),
                         paste0(results_folder_base,"CT4"))




plot_folder_base <- "plots/TIGER_tissue_100_GT"

for(max_GT_depth in 1:5){
  plot_folder_base_GT <- paste0(plot_folder_base,max_GT_depth,"/")
  
  ct_plot_folder = list(paste0(plot_folder_base_GT,"all/"),
                        paste0(plot_folder_base_GT,"CT1/"),
                        paste0(plot_folder_base_GT,"CT2/"),
                        paste0(plot_folder_base_GT,"CT3/"),
                        paste0(plot_folder_base_GT,"CT4/"))
  
  
  
  lapply(ct_plot_folder,function(f)dir.create(f,recursive = TRUE))
  
  
  
  for(i_analysis in 1:length(ct_results_folder)){
    
    results_folder = ct_results_folder[[i_analysis]]
    plot_folder = ct_plot_folder[[i_analysis]]
    
    current_cell_type = c("all","CT1","CT2","CT3","CT4")[[i_analysis]]
    
    # Analysis
    results <- list.dirs(results_folder, recursive = FALSE)
    misty.results <- collect_results(results)
    
    
    # Plot contributions and improvement
    misty.results %>% plot_view_contributions() 
    ggsave(paste0(plot_folder,"view_contribution.pdf"), width = 10,height = 3.2)
    
    # difference between multi-intra 
    misty.results %>% plot_improvement_stats()
    ggsave(paste0(plot_folder,"r2_improvement.pdf"), width = 10,height = 3.2)
    
    
    # variance explained by intra-view only. 
    misty.results %>% plot_improvement_stats(measure = "intra.R2")
    ggsave(paste0(plot_folder,"r2_by_intra.pdf"), width = 10,height = 3.2)
    
    
    
    #misty.results %>% plot_improvement_stats(measure = "multi.R2")
    write_tsv(misty.results$improvements.stats, paste0(plot_folder,"improvements.stats.txt"))
    
    
    # Plot aggregated importances
    misty.results %>% plot_interaction_heatmap(view = "intra",cutoff =  0.5)
    ggsave(paste0(plot_folder,"interaction_heatmap_intra.pdf"), width = 10,height = 10)
    
    misty.results %>% plot_interaction_heatmap("para.10", 0.5) 
    ggsave(paste0(plot_folder,"interaction_heatmap_para.pdf"), width = 10,height = 10)
    
    misty.results %>% plot_contrast_heatmap("intra", "para.10", 0.5)
    ggsave(paste0(plot_folder,"contrast.pdf"), width = 10,height = 10)
    
    #misty.results %>% plot_interaction_communities(view = "intra")
    #needs manual save
    
    #misty.results %>% plot_interaction_communities(view = "para.10", cutoff = 0.5)
    #needs manual save
    # read_csv("data/synthetic/true_connections_nolig.csv")
    
    
    
    
    true.interactions <- read_csv("data/TIGER_tissue/cell_type_interactions.csv") %>%
      mutate(Vmax = ifelse(from==target,NA,Vmax))
    all_markers = misty.results$improvements$target %>% unique()
    
    # get the ground truth
    # starts from the interactions 
    true.connections_directed <- get_connections(current_cell_type, true.interactions, all_markers)
    
    true.connections_directed <- true.connections_directed %>% filter(view=="intra" | is.na(depth) | (depth < max_GT_depth))
    
    true.connections_reversed <- true.connections_directed %>%
      filter(present) %>%
      filter(node1!=node2) %>%
      mutate(nodeX = node1,
             node1 = node2,
             node2 = nodeX) %>%
      select(-nodeX)
    
    true.connections <- bind_rows(true.connections_directed,true.connections_reversed) %>%
      group_by(node1,node2,view,cell_type) %>%
      summarise(
        present = any(present)) %>% ungroup()
    
    
    # Plot true connections
    ggtrue.intra <- ggplot(true.connections %>% filter(view == "intra")) + 
      geom_tile(aes(x = node1, y = node2, fill = as.factor(present))) +
      scale_fill_discrete(type = c("white", muted("blue"))) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) + ggtitle("true.intra")
    
    ggtrue.para <- ggplot(true.connections %>% filter(view == "para")) + 
      geom_tile(aes(x = node1, y = node2, fill = as.factor(present))) +
      scale_fill_discrete(type =  c("white", muted("blue"))) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) + ggtitle("true.para")
    
    
    plot_grid(ggtrue.intra, ggtrue.para)
    ggsave(paste0(plot_folder,"gold_standard.pdf"),width = 20,height = 10)
    
    
    # Distributions of AUROC and AUPRC
    
    
    joined_intra <- true.connections %>% filter(view == "intra") %>% 
      left_join(misty.results$importances %>% filter(view == "intra") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    
    intra.roc <- joined_intra %>% group_by(sample) %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    intra.pr <- joined_intra %>% group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    joined_para <- true.connections %>% filter(view == "para") %>% 
      left_join(misty.results$importances %>% filter(view == "para.10") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    
    para.roc <- joined_para %>% group_by(sample) %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    para.pr <- joined_para %>% group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    # para with pre-filtering: 
    trim = 0
    targets <- misty.results$improvements.stats %>% filter(measure == "gain.R2", mean > trim) %>% pull(target)
    
    gain <- misty.results$improvements.stats %>% filter(measure == "gain.R2")
    
    joind_para_trim =joined_para %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE)))  
    
    para.roc_trim <- joind_para_trim  %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    para.pr_trim <- joined_para %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE))) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    
    
    viodata.roc <- data.frame(intra.roc = intra.roc %>% unlist,
                              para.roc = para.roc %>% unlist,
                              para.roc_trim = para.roc_trim ) %>% 
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.roc, aes(x=Type, y=AUC)) + geom_violin(aes(fill=Type)) + 
      geom_jitter() + 
      geom_hline(yintercept=0.5, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=0.5, color="#00BFC4", linetype="dashed") +
      theme(axis.text = element_text(family = "Arial")) +
      theme_classic() + ylim(0.0,1)
    
    ggsave(paste0(plot_folder,"AUC_ROC_intra_para.pdf"), width = 4,height = 3.2)
    
    intra.intercept <- sum(true.connections %>% filter(view == "intra", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "intra", !is.na(present)) %>% nrow())
    para.intercept <- sum(true.connections %>% filter(view == "para", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "para", !is.na(present)) %>% nrow())
    
    viodata.pr <- data.frame(intra.pr = intra.pr %>% unlist,
                             para.pr = para.pr %>% unlist,
                             para.pr_trim = para.pr_trim) %>%
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.pr, aes(x=Type, y=AUC)) + geom_violin(aes(fill=Type)) + 
      geom_jitter() + 
      geom_hline(yintercept=intra.intercept, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=para.intercept, color="#00BFC4", linetype="dashed") +
      theme_classic() + ylim(0,1)
    
    ggsave(paste0(plot_folder,"AUC_PR_intra_para.pdf"), width = 4,height = 3.2)
    
    
    # Plot aggregated ROC and PR curves
    
    # for iso-F1 and iso-J lines
    r <- seq(0,1,by=1e-3)
    ks <- c(0.1,0.2,0.5,0.8)
    para.col = "#EC008C"
    intra.col = "#00A651"
    para_trim.col = "#377EB8"
    
    pdf(file = paste0(plot_folder,"ROC_para_intra.pdf"),width = 10,height = 10)
    tidy.intra.agg <- misty.results$importances.aggregated %>% filter(view == "intra")
    joined.intra.agg <- true.connections %>% filter(view == "intra") %>% 
      left_join(tidy.intra.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_intra <- roc.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                           curve = TRUE, rand.compute = TRUE)
    plot(roc_intra, color=intra.col, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE)
    
    # add para
    tidy.para.agg <- misty.results$importances.aggregated %>% filter(view == "para.10")
    joined.para.agg <- true.connections %>% filter(view == "para") %>% 
      left_join(tidy.para.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_para = roc.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                         curve = TRUE, rand.compute = TRUE)
    plot(roc_para, color=para.col, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    # add para trim
    
    tidy.para_trim.agg <- misty.results$importances.aggregated %>% filter(view == "para.10") %>%
      mutate(Importance = ifelse(Target %in% targets, Importance, min(Importance,na.rm = TRUE)))
    
    joined.para_trim.agg <- true.connections %>% filter(view == "para") %>% 
      left_join(tidy.para_trim.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_para_trim = roc.curve(joined.para_trim.agg %>% pull(Importance), weights.class0 = joined.para_trim.agg %>% pull(present), 
                              curve = TRUE, rand.compute = TRUE)
    plot(roc_para_trim, color=para_trim.col, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    
    
    text(.4,.1,labels = paste("AUC(para)=", round(roc_para$auc,digits = 3)),adj = 0,col=para.col)
    text(.4,.2,labels = paste("AUC(intra)=", round(roc_intra$auc,digits = 3)),adj = 0,col=intra.col)
    text(.4,.3,labels = paste("AUC(para_trim)=", round(roc_para_trim$auc,digits = 3)),adj = 0,col=para_trim.col)
    
    ks %>% walk(function(k){
      s <- r + k
      s[s>1] <- NA
      lines(r, s, ylim=c(0,1),xlim=c(0,1), col="gray80")
    })
    dev.off()
    
    
    
    # PR curves
    para.col = "#EC008C"
    intra.col = "#00A651"
    
    pdf(file = paste0(plot_folder,"PR_para_intra.pdf"),width = 10,height = 10)
    pr_intra <- pr.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                         curve = TRUE)
    
    plot(pr_intra, color=intra.col, auc.main=FALSE)
    
    
    pr_para <- pr.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                        curve = TRUE)
    plot(pr_para, color=para.col, rand.plot = TRUE, auc.main=FALSE, add = TRUE)
    
    
    pr_para_trim <- pr.curve(joined.para_trim.agg %>% pull(Importance), weights.class0 = joined.para_trim.agg %>% pull(present), 
                             curve = TRUE)
    plot(pr_para_trim, color=para_trim.col, rand.plot = TRUE, auc.main=FALSE, add = TRUE)
    
    
    lines(x = c(0,1), y = c(intra.intercept,intra.intercept), lty = 2, col = intra.col)
    lines(x = c(0,1), y = c(para.intercept, para.intercept), lty = 2, col = para.col)
    lines(x = c(0,1), y = c(para.intercept, para.intercept), lty = 2, col = para.col)
    
    ks %>% walk(function(k){
      p <- k*r/(2*r - k)
      p[p<=0 | p > 1.5] <- NA
      lines(r, p,ylim=c(0,1),xlim=c(0,1), col="gray80")
    })
    
    text(.4,.9,labels = paste("AUC(para_trim)=", round(pr_para_trim$auc.integral,digits = 3)),adj = 0,col=para_trim.col)
    text(.4,.8,labels = paste("AUC(para)=", round(pr_para$auc.integral,digits = 3)),adj = 0,col=para.col)
    text(.4,.7,labels = paste("AUC(intra)=", round(pr_intra$auc.integral,digits = 3)),adj = 0,col=intra.col)
    
    
    dev.off()
    
  }
}



####### ADD Recall - sensitivity curves


ct_results_folder = list(paste0(results_folder_base,"all"),
                         paste0(results_folder_base,"CT1"),
                         paste0(results_folder_base,"CT2"),
                         paste0(results_folder_base,"CT3"),
                         paste0(results_folder_base,"CT4"))




plot_folder_base <- "plots/TIGER_tissue_100_GT"
i_list = 0
PRlist = list()


for(i_analysis in 1:length(ct_results_folder)){
  
  results_folder = ct_results_folder[[i_analysis]]
  
  
  current_cell_type = c("all","CT1","CT2","CT3","CT4")[[i_analysis]]
  
  # Analysis
  results <- list.dirs(results_folder, recursive = FALSE)
  misty.results <- collect_results(results)
  
  
  true.interactions <- read_csv("data/TIGER_tissue/cell_type_interactions.csv") %>%
    mutate(Vmax = ifelse(from==target,NA,Vmax))
  all_markers = misty.results$improvements$target %>% unique()
  
  
  for(max_GT_depth in 1:5){
    
    
    # get the ground truth
    # starts from the interactions 
    true.connections_directed <- get_connections(current_cell_type, true.interactions, all_markers)
    
    true.connections_directed <- true.connections_directed %>% filter(view=="intra" | is.na(depth) | (depth < max_GT_depth))
    
    
    
    bi_directional = FALSE
    if(bi_directional){
      true.connections_reversed <- true.connections_directed %>%
        filter(present) %>%
        filter(node1!=node2) %>%
        mutate(nodeX = node1,
               node1 = node2,
               node2 = nodeX) %>%
        select(-nodeX)
      
      true.connections <- bind_rows(true.connections_directed,true.connections_reversed) %>%
        group_by(node1,node2,view,cell_type) %>%
        summarise(
          present = any(present)) %>% ungroup()
      
    }else{
      true.connections <- true.connections_directed
    }
    
    
    # Distributions of AUROC and AUPRC
    joined_para <- true.connections %>% filter(view == "para") %>% 
      left_join(misty.results$importances %>% filter(view == "para.10") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    
    para.pr <- joined_para %>% group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral,
                PR  =  list(pr.curve(Importance, weights.class0 = present, curve = TRUE)))
    
    pr_curve <-  para.pr %>% mutate(pr_curve = map(PR,function(x)as_tibble(x$curve))) %>% unnest(pr_curve) %>%
      select(-PR,-auc)
    
    colnames(pr_curve) = c("sample","Recall","Precision","threshold")
    
    
    
    # para with pre-filtering: 
    trim = 0
    targets <- misty.results$improvements.stats %>% filter(measure == "gain.R2", mean > trim) %>% pull(target)
    gain <- misty.results$improvements.stats %>% filter(measure == "gain.R2")
    
    
    para.pr_trim <- joined_para %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE))) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral,
                PR  =  list(pr.curve(Importance, weights.class0 = present, curve = TRUE)))
    
    
    pr_curve_trim <-  para.pr_trim %>% mutate(pr_curve = map(PR,function(x)as_tibble(x$curve))) %>% unnest(pr_curve) %>%
      select(-PR,-auc)
    
    colnames(pr_curve_trim) = c("sample","Recall","Precision","threshold")
    
    # random model: 
    
    joined_para <- true.connections %>% filter(view == "para") %>% 
      left_join(misty.results$importances %>% filter(view == "para.10") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    # shuffling the imporntance
    set.seed(156)
    para.pr_random <- joined_para %>% group_by(sample,node2) %>%
      mutate(Importance = sample(Importance)) %>%
      group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral,
                PR  =  list(pr.curve(Importance, weights.class0 = present, curve = TRUE)))
    
    pr_curve_random <-  para.pr_random %>% mutate(pr_curve = map(PR,function(x)as_tibble(x$curve))) %>% unnest(pr_curve) %>%
      select(-PR,-auc)
    
    colnames(pr_curve_random) = c("sample","Recall","Precision","threshold")
    
    # random normal model: 
    
    # sampling informance from standard normal
    set.seed(156)
    para.pr_randomSN <- joined_para %>% group_by(sample,node2) %>%
      mutate(Importance = rnorm(length(Importance))) %>%
      group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral,
                PR  =  list(pr.curve(Importance, weights.class0 = present, curve = TRUE)))
    
    pr_curve_randomSN <-  para.pr_randomSN %>% mutate(pr_curve = map(PR,function(x)as_tibble(x$curve))) %>% unnest(pr_curve) %>%
      select(-PR,-auc)
    
    colnames(pr_curve_randomSN) = c("sample","Recall","Precision","threshold")
    
    
    
    # Put things together
    
    pr_curve$method = "no_trim"
    pr_curve_trim$method = "trim0"
    pr_curve_random$method = "random"
    pr_curve_randomSN$method = "SNrandom"
    
    i_list = i_list + 1 
    PRlist[[i_list]] <- bind_rows(pr_curve,pr_curve_trim,pr_curve_random,pr_curve_randomSN) %>%
      mutate(cell_type = current_cell_type,
             GT_depth = max_GT_depth)
  }
}

PRall <- PRlist %>% bind_rows()
PRall %>% ggplot(aes(threshold,Recall)) + geom_line(aes(col = sample)) + facet_grid(GT_depth~method) + guides(col = "none")


PRall %>%
  filter(grepl("IMC3",sample)) %>%
  ggplot(aes(threshold,Recall)) + geom_line(aes(col = as.factor(method))) + facet_grid(cell_type~GT_depth) 


# strongest false positives:


joined_para %>% filter(present == FALSE) %>% arrange(desc(Importance)) 



raw_imp <- read_csv("./results/TIGER_tissue_100/CT4/IMC2/importances_X103_para.10.txt")

raw_imp %>% arrange(desc(imp))
# # A tibble: 96 × 2
# target       imp
# <chr>      <dbl>
#   1 R5     0.000153 
# 2 X63    0.000130 
# 3 X58    0.000123 

raw_imp$imp %>% mean
# [1] 7.394136e-05



raw_imp <- read_csv("./results/TIGER_tissue_100/CT4/IMC2/importances_R5_para.10.txt")

raw_imp %>% arrange(desc(imp))
# A tibble: 96 × 2
# target   imp
# <chr>  <dbl>
#   1 X134    36.4
# 2 X29     35.2
# 3 X91     34.3

raw_imp$imp %>% hist()

raw_imp$imp %>% mean
# [1]24.56677




# meanwhile something that actuially has effect: 
# X47 is direct parent of L3
raw_impX47 <- read_csv("./results/TIGER_tissue_100/CT4/IMC2/importances_X47_para.10.txt")

raw_impX47 %>% arrange(desc(imp))
# A tibble: 96 × 2
# target   imp
# <chr>  <dbl>
#   1 X91     4.01
# 2 X25     3.72
# 3 R5      3.66
# 4 X124    3.64
raw_impX47$imp %>% hist()


flist <- list.files("./results/TIGER_tissue_100/CT2/IMC2/",pattern = "para",full.names = TRUE)
feature <- strsplit(flist,split = "_") %>% lapply(function(x)x[4]) %>% unlist()
importances <- lapply(flist,read_csv)
names(importances) <- feature
importances <- bind_rows(importances,.id ="feature")


importances %>% filter(target == "X103") %>% arrange(desc(imp))

