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


## This is the script producing the figures in the revision for Genome Bioology.

input_folder = "./data/TIGER_general_network_v8/"
results_folder_base = "results/TIGER_general_network_v8/"
plot_folder_base = "plots/TIGER_general_network_v8_revision"

# Run MISTy
data <- list.dirs(input_folder, recursive = FALSE)[3:4]
plan(multisession)

l <- 10

exp_id = "small_network"

# I first eliminate those targets which seem to be not expressed in the cell type
d = data[[1]]
measurements <- read_csv(paste0(d, "/small_network_position_expression_real_cells.csv")) %>%
  select(-starts_with("L"))

nodes <- colnames(measurements)[3:(ncol(measurements)-1)]
set.seed(789)


observed_nodes <- nodes
  # nodes[!nodes %>% grepl("^R",x = .)] %>% sample(.,10)
  # setdiff(nodes,c("X10","X22","X4","X8","X16","X20","X21","X8","X17",
  #                                "X25","X31","X23","X29","X28"))

data %>% walk(function(d) {
  print(paste("processing:",d))
  
  all <- read_csv(paste0(d, "/small_network_position_expression_real_cells.csv")) %>%
    select(-starts_with("L")) %>%
    select(row,col,type,all_of(observed_nodes))
  
  
  expr <- all %>% select(-row, -col,-type)

  pos <- all %>% select(row, col)
  
  views <- create_initial_view(expr) %>% select_markers("intraview", where(~sd(.) > .1)) %>%  add_paraview(pos, l)
  
  run_misty(views, results.folder = paste0(results_folder_base,
                                           "all/",
                                           str_extract(d, "IMC[0-9]+"), "/"
  ), #bypass.intra = TRUE
  )
  
  types <- all %>% pull(type)
  
  unique(types) %>% walk(function(t) {
    
    # - filtering the intraview for expressed markers (sd > 0.1)
    # - these markers are also used as targets in the misty model by default,
    # so this way misty does not build models for non-expressed markers
    filtered_ct_view <- views %>% filter_views(which(types == t)) %>%
      select_markers("intraview", where(~sd(.) > .1)) 
    
                                      
    run_misty(filtered_ct_view,
              results.folder = paste0(
                results_folder_base, t, "/",
                str_extract(d, "IMC[0-9]+"), "/"
              ), #bypass.intra = TRUE
    )
  })
})




ct_results_folder = list(paste0(results_folder_base,"all"),
                         paste0(results_folder_base,"CT1"),
                         paste0(results_folder_base,"CT2"),
                         paste0(results_folder_base,"CT3"),
                         paste0(results_folder_base,"CT4"))






for(max_GT_depth in c(1)){
  plot_folder_base_GT <- paste0(plot_folder_base,"_revision_GT",max_GT_depth,"/")
  
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
    results <- list.dirs(results_folder, recursive = FALSE)[3:4]
    misty.results <- collect_results(results)
    
    
    # Plot contributions and improvement
    misty.results %>% plot_view_contributions()
    ggsave(paste0(plot_folder,"view_contribution.pdf"), width = 5,height = 3.2)
    
    # difference between multi-intra 
    misty.results %>% plot_improvement_stats()
    ggsave(paste0(plot_folder,"r2_improvement.pdf"), width = 5,height = 3.2)
    
    
    # variance explained by intra-view only. 
    misty.results %>% plot_improvement_stats(measure = "intra.R2")
    ggsave(paste0(plot_folder,"r2_by_intra.pdf"), width = 5,height = 3.2)
    
    #misty.results %>% plot_improvement_stats(measure = "multi.R2")
    write_tsv(misty.results$improvements.stats, paste0(plot_folder,"improvements.stats.txt"))
    
    # Plot aggregated importances
    misty.results %>% plot_interaction_heatmap(view = "intra",cutoff =  0.5)
    ggsave(paste0(plot_folder,"interaction_heatmap_intra.pdf"), width = 5,height = 5)
    
    misty.results %>% plot_interaction_heatmap("para.10", 0.5) 
    ggsave(paste0(plot_folder,"interaction_heatmap_para.pdf"), width = 5,height = 5)
    
    
    true.interactions <- read_csv(file.path(input_folder,"small_network_cell_type_interactions.csv")) %>%
      mutate(Vmax = ifelse(from==target,NA,Vmax))
    all_markers = misty.results$improvements$target %>% unique()
    
    # get the ground truth
    # starts from the interactions 
    true.connections_directed <- get_connections(current_cell_type, true.interactions, all_markers)
    true.connections_directed <- true.connections_directed %>% 
      mutate(present = ifelse(present == TRUE & depth < max_GT_depth, TRUE, FALSE)) 
    # %>% 
    #   filter(view=="intra" | is.na(depth) | (depth < max_GT_depth))
    
    bidirectional = TRUE
    if(bidirectional){
      
      true.connections_reversed <- true.connections_directed %>%
        filter(present) %>%
        filter(node1!=node2) %>%
        mutate(nodeX = node1,
               node1 = node2,
               node2 = nodeX) %>%
        select(-nodeX)
      
      true.connections_all <- bind_rows(true.connections_directed,true.connections_reversed) %>%
        group_by(node1,node2,view,cell_type) %>%
        summarise(
          present = any(present)) %>% ungroup()
    }else{
      true.connections_all <- true.connections_directed
    }
    # keep only the interactions which are observed
    true.connections <- true.connections_all %>%
      filter(node1 %in% observed_nodes, node2 %in% observed_nodes) %>%
      mutate(present = ifelse(node1==node2,NA,present))
    
    # Plot true connections
    ggtrue.intra <- ggplot(true.connections %>% filter(view == "intra")) + 
      geom_tile(aes(x = node1, y = node2, fill = as.factor(present))) +
      scale_fill_discrete(type = c("white", muted("blue"))) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank()) +
      ggtitle("true.intra") + guides(fill="none")+
      xlab("Predictor")+ ylab("Target")
    
    ggtrue.para <- ggplot(true.connections %>% filter(view == "para")) + 
      geom_tile(aes(x = node1, y = node2, fill = as.factor(present))) +
      scale_fill_discrete(type =  c("white", muted("blue"))) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90), 
            panel.grid = element_blank()) + ggtitle("true.para")+ guides(fill="none")+
    xlab("Predictor") + ylab("Target")
    
    plot_grid(ggtrue.intra, ggtrue.para)
    ggsave(paste0(plot_folder,"gold_standard.pdf"),width = 5,height = 5)
    
    
    # Distributions of AUROC and AUPRC
    
    
    joined_intra <- true.connections %>% filter(view == "intra") %>% 
      left_join(misty.results$importances %>% filter(view == "intra"), by = c("node1" = "Predictor", "node2" = "Target")) %>% 
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


    # para.roc <- joined_para %>% group_by(sample) %>%
    #   summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
    #   pull(auc)
    # 
    # para.pr <- joined_para %>% group_by(sample) %>%
    #   summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
    #   pull(auc)
    
    # para with pre-filtering: 
    trim = 0
    targets <- misty.results$improvements.stats %>% filter(measure == "gain.R2", mean > trim) %>% pull(target)
    
    gain <- misty.results$improvements.stats %>% filter(measure == "gain.R2")
    
    joined_para =joined_para %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE)))  
    
    para.roc<- joined_para  %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    para.pr <- joined_para %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE))) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    
    
    viodata.roc <- data.frame(intra.roc = intra.roc %>% unlist,
                              para.roc = para.roc %>% unlist) %>% 
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.roc, aes(x=Type, y=AUC)) + geom_boxplot() + 
      geom_jitter(aes(col=Type),width = 0.1) + 
      geom_hline(yintercept=0.5, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=0.5, color="#00BFC4", linetype="dashed") +
      theme_classic() +
      ylim(0.0,1) + ggtitle("Area under ROC curve") + guides(col="none") + xlab("")
    
    ggsave(paste0(plot_folder,"AUC_ROC_intra_para.pdf"), width = 2.5,height = 3)
    
    intra.intercept <- sum(true.connections %>% filter(view == "intra", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "intra", !is.na(present)) %>% nrow())
    para.intercept <- sum(true.connections %>% filter(view == "para", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "para", !is.na(present)) %>% nrow())
    
    viodata.pr <- data.frame(intra.pr = intra.pr %>% unlist,
                             para.pr = para.pr %>% unlist
                             ) %>%
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.pr, aes(x=Type, y=AUC)) + geom_boxplot() + 
      geom_jitter(aes(col=Type),width = 0.1) + 
      geom_hline(yintercept=intra.intercept, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=para.intercept, color="#00BFC4", linetype="dashed") +
      theme_classic() + ylim(0,1) +ggtitle("Area under PR curve") + guides(col="none") + xlab("")
    
    ggsave(paste0(plot_folder,"AUC_PR_intra_para.pdf"), width = 2.5,height = 3)
    
    
    # Plot aggregated ROC and PR curves
    
    # for iso-F1 and iso-J lines
    r <- seq(0,1,by=1e-3)
    ks <- c(0.1,0.2,0.5,0.8)
    para.col = "#EC008C"
    intra.col = "#00A651"
   
    
    pdf(file = paste0(plot_folder,"ROC_para_intra.pdf"),width = 3,height = 3)
   
     tidy.intra.agg <- misty.results$importances.aggregated %>% filter(view == "intra")
    joined.intra.agg <- true.connections %>% filter(view == "intra") %>% 
      left_join(tidy.intra.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_intra <- roc.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                           curve = TRUE, rand.compute = TRUE)
    plot(roc_intra, color=intra.col, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE)
    
    # add para
    tidy.para.agg <- misty.results$importances.aggregated %>% filter(view == "para.10") %>%
      mutate(Importance = ifelse(Target %in% targets, Importance, min(Importance,na.rm = TRUE)))
    joined.para.agg <- true.connections %>% filter(view == "para") %>% 
      left_join(tidy.para.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_para = roc.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                         curve = TRUE, rand.compute = TRUE)
    plot(roc_para, color=para.col, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    
    text(.2,.1,labels = paste("AUC(para)=", round(roc_para$auc,digits = 3)),adj = 0,col=para.col)
    text(.2,.25,labels = paste("AUC(intra)=", round(roc_intra$auc,digits = 3)),adj = 0,col=intra.col)
    
    
    ks %>% walk(function(k){
      s <- r + k
      s[s>1] <- NA
      lines(r, s, ylim=c(0,1),xlim=c(0,1), col="gray80")
    })
    dev.off()
    
    
    
    # PR curves
    para.col = "#EC008C"
    intra.col = "#00A651"
    
    pdf(file = paste0(plot_folder,"PR_para_intra.pdf"),width = 3,height = 3)
    pr_intra <- pr.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                         curve = TRUE)
    
    plot(pr_intra, color=intra.col, auc.main=FALSE)
    
    
    pr_para <- pr.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                        curve = TRUE)
    plot(pr_para, color=para.col, rand.plot = TRUE, auc.main=FALSE, add = TRUE)
    
    lines(x = c(0,1), y = c(intra.intercept,intra.intercept), lty = 2, col = intra.col)
    lines(x = c(0,1), y = c(para.intercept, para.intercept), lty = 2, col = para.col)
    
    
    ks %>% walk(function(k){
      p <- k*r/(2*r - k)
      p[p<=0 | p > 1.5] <- NA
      lines(r, p,ylim=c(0,1),xlim=c(0,1), col="gray80")
    })
    
    text(.2,.8,labels = paste("AUC(para)=", round(pr_para$auc.integral,digits = 3)),adj = 0,col=para.col)
    text(.2,.6,labels = paste("AUC(intra)=", round(pr_intra$auc.integral,digits = 3)),adj = 0,col=intra.col)
    
    
    dev.off()
    
  }
}





# general and Celltype 1 network together
max_GT_depth = 1

ct_plot_folder = "plots/TIGER_general_network_V8_GT_revision_GT1_allCT1/"


dir.create(ct_plot_folder,recursive = TRUE)

results_folder_all = ct_results_folder[[1]]
results_folder_CT1 = ct_results_folder[[2]]
plot_folder = ct_plot_folder

current_cell_type = c("all","CT1")

# Analysis
results <- list.dirs(results_folder_all, recursive = FALSE)[3:4]
misty.results_all <- collect_results(results)

results <- list.dirs(results_folder_CT1, recursive = FALSE)[3:4]
misty.results_CT1 <- collect_results(results)


    # Plot contributions and improvement
bind_rows(misty.results_all$contributions.stats %>% add_column(case = "no CT"), 
          misty.results_CT1$contributions.stats %>% add_column(case = "CT1") )%>% 
  ggplot2::ggplot(aes(x = target,  y = fraction)) + 
  ggplot2::geom_col(aes(group = view, fill = view)) + 
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::theme_classic() + ggplot2::ylab("Contribution") + 
  ggplot2::xlab("Target") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,  hjust = 1,size = 7)) +
  facet_wrap(~case, scales = "free_x")

    ggsave(paste0(plot_folder,"view_contribution.pdf"), width = 6,height =2.8)
    
    # difference between multi-intra 
    # Plot contributions and improvement
    set2.orange <- "#FC8D62"
    set2.purple <- "#8DA0CB"
    tmp = bind_rows(misty.results_all$improvements.stats  %>% add_column(case = "no CT"), 
              misty.results_CT1$improvements.stats  %>% add_column(case = "CT1") ) %>%
      filter(measure == "gain.R2") 
    pos_r2_targets <- tmp %>% group_by(target)  %>% summarise(max_r2 = max(mean)) %>% filter(max_r2>0.01) %>% pull(target)
    
    tmp %>%
      filter(target %in% pos_r2_targets) %>% 
      ggplot2::ggplot(aes(x = stats::reorder(.data$target, 
                                             -.data$mean), y = .data$mean)) + 
      ggplot2::geom_pointrange(aes(ymin = .data$mean - .data$sd, ymax = .data$mean + .data$sd)) +
      ggplot2::geom_point(aes(color = case)) + 
      ggplot2::theme_classic() + ggplot2::ylab("gain.R2") + ggplot2::xlab("Target") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                         hjust = 1,size = 9)) 
    
    ggsave(paste0(plot_folder,"r2_improvement.pdf"), width = 4,height = 2.8)
    
    
    # variance explained by intra-view only. 
    bind_rows(misty.results_all$improvements.stats  %>% add_column(case = "no CT"), 
              misty.results_CT1$improvements.stats  %>% add_column(case = "CT1") )%>% 
      filter(measure == "intra.R2") %>%
      filter(mean > 0) %>% 
      ggplot2::ggplot(aes(x = stats::reorder(.data$target, 
                                             -.data$mean), y = .data$mean)) + 
      ggplot2::geom_pointrange(aes(ymin = .data$mean - .data$sd, ymax = .data$mean + .data$sd)) +
      ggplot2::geom_point(aes(color = case)) + 
      ggplot2::theme_classic() + ggplot2::ylab("intra.R2") + ggplot2::xlab("Target") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                         hjust = 1,size=7)) 
    
    ggsave(paste0(plot_folder,"r2_by_intra.pdf"), width = 4,height = 2.2)
    
    
    # Plot aggregated importances
    misty.results_all %>% plot_interaction_heatmap(view = "intra",cutoff =  0.5)
    ggsave(paste0(plot_folder,"interaction_heatmap_intra_all.pdf"), width = 5,height = 5)
    
    misty.results_CT1 %>% plot_interaction_heatmap(view = "intra",cutoff =  0.5)
    ggsave(paste0(plot_folder,"interaction_heatmap_intra_CT1.pdf"), width = 5,height = 5)
    
    
    misty.results_all %>% plot_interaction_heatmap("para.10", 0.5) 
    ggsave(paste0(plot_folder,"interaction_heatmap_para_all.pdf"), width = 5,height = 5)
    
    misty.results_CT1 %>% plot_interaction_heatmap("para.10", 0.5) 
    ggsave(paste0(plot_folder,"interaction_heatmap_para_CT1.pdf"), width = 5,height = 5)
    
    
    true.interactions <- read_csv(file.path(input_folder,"small_network_cell_type_interactions.csv")) %>%
      mutate(Vmax = ifelse(from==target,NA,Vmax))
    all_markers = misty.results$improvements$target %>% unique()
    
    # get the ground truth
    # starts from the interactions 
    wrap_GT <- function(cct,true.interactions,all_markers,max_GT_depth){
      true.connections_directed <- get_connections(cct, true.interactions, all_markers)
      true.connections_directed <- true.connections_directed %>% 
        mutate(present = ifelse(present == TRUE & depth < max_GT_depth, TRUE, FALSE)) 
      # %>% 
      #   filter(view=="intra" | is.na(depth) | (depth < max_GT_depth))
      
      bidirectional = TRUE
      if(bidirectional){
        
        true.connections_reversed <- true.connections_directed %>%
          filter(present) %>%
          filter(node1!=node2) %>%
          mutate(nodeX = node1,
                 node1 = node2,
                 node2 = nodeX) %>%
          select(-nodeX)
        
        true.connections_all <- bind_rows(true.connections_directed,true.connections_reversed) %>%
          group_by(node1,node2,view,cell_type) %>%
          summarise(
            present = any(present)) %>% ungroup()
      }else{
        true.connections_all <- true.connections_directed
      }
      # keep only the interactions which are observed
      true.connections <- true.connections_all %>% filter(node1 %in% observed_nodes, node2 %in% observed_nodes)%>%
        mutate(present = ifelse(node1==node2,NA,present))
      return(true.connections)
    }
    
    true.connections_all <- wrap_GT(cct = "all",true.interactions,all_markers,max_GT_depth)
    true.connections_CT1 <- wrap_GT(cct = "CT1",true.interactions,all_markers,max_GT_depth)
    
    
    
    set2.blue <- "#8DA0CB"
    # Plot true connections
    ggtrue.intra <- misty.results_all$importances.aggregated %>%
      dplyr::filter(view=="intra") %>%
      ggplot( ggplot2::aes(x = .data$Predictor, 
                                                   y = .data$Target)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) + 
      ggplot2::scale_fill_gradient2(low = "white", mid = "white", 
                                    high = set2.blue, midpoint = 0.5) +
      geom_point(data = true.connections_all %>% 
                   filter(view == "intra",
                          present==TRUE),
                 aes(x = node1, y = node2),col= set2.orange,shape =4,size = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank()) +
      ggtitle("true.intra") + guides(fill="none")+
      xlab("Predictor")+ ylab("Target")
    
    ggtrue.para <- misty.results_all$importances.aggregated %>%
      dplyr::filter(view=="para.10") %>%
      ggplot( ggplot2::aes(x = .data$Predictor, 
                           y = .data$Target)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) + 
      ggplot2::scale_fill_gradient2(low = "white", mid = "white", 
                                    high = set2.blue, midpoint = 0) +
      geom_point(data = true.connections_all %>% 
                   filter(view == "para",
                          present==TRUE),
                 aes(x = node1, y = node2),col= set2.orange,shape =4,size = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank()) +
      ggtitle("true.para") +
      xlab("Predictor")+ ylab("Target")
    
    plot_grid(ggtrue.intra, ggtrue.para)
    ggsave(paste0(plot_folder,"prediction_with_gold_standard_all_legend.pdf"),width = 8,height = 4)
    #ggsave(paste0(plot_folder,"prediction_with_gold_standard_all.pdf"),width = 8,height = 4)
    
    
    
    set2.blue <- "#8DA0CB"
    # Plot true connections
    ggtrue.intra <- misty.results_CT1$importances.aggregated %>%
      dplyr::filter(view=="intra") %>%
      ggplot( ggplot2::aes(x = .data$Predictor, 
                           y = .data$Target)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) + 
      ggplot2::scale_fill_gradient2(low = "white", mid = "white", 
                                    high = set2.blue, midpoint = 0.5) +
      geom_point(data = true.connections_CT1 %>% 
                   filter(view == "intra",
                          present==TRUE),
                 aes(x = node1, y = node2),col= set2.orange,shape =4,size = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank()) +
      ggtitle("true.intra") + guides(fill="none")+
      xlab("Predictor")+ ylab("Target")
    
    ggtrue.para <- misty.results_CT1$importances.aggregated %>%
      dplyr::filter(view=="para.10") %>%
      ggplot( ggplot2::aes(x = .data$Predictor, 
                           y = .data$Target)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance)) + 
      ggplot2::scale_fill_gradient2(low = "white", mid = "white", 
                                    high = set2.blue, midpoint = 0) +
    geom_point(data = true.connections_CT1 %>% 
                   filter(view == "para",
                          present==TRUE),
                 aes(x = node1, y = node2),col= set2.orange,shape =4,size = 4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank()) +
      ggtitle("true.para") + guides(fill="none")+
      xlab("Predictor")+ ylab("Target")
    
    plot_grid(ggtrue.intra, ggtrue.para)
    ggsave(paste0(plot_folder,"prediction_with_gold_standard_CT1.pdf"),width = 8,height = 4)
    
    
    # Distributions of AUROC and AUPRC
    
    # all
    joined_intra_all <- true.connections_all %>% filter(view == "intra") %>% 
      left_join(misty.results_all$importances %>% filter(view == "intra") %>% mutate(method = "no CT") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    intra.roc_all <- joined_intra_all %>% group_by(sample) %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    intra.pr_all <- joined_intra_all %>% group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)

    
    # CT1
    joined_intra_CT1 <- true.connections_CT1 %>% filter(view == "intra") %>% 
      left_join(misty.results_CT1$importances %>% filter(view == "intra") %>% mutate(method = "no CT") , by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    
    intra.roc_CT1 <- joined_intra_CT1 %>% group_by(sample) %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    intra.pr_CT1 <- joined_intra_CT1 %>% group_by(sample) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    

    
    # para with pre-filtering ALL 
    trim = 0
    targets <- misty.results_all$improvements.stats %>% filter(measure == "gain.R2", mean > trim) %>% pull(target)
    
    gain <- misty.results_all$improvements.stats %>% filter(measure == "gain.R2")
    
    joined_para_all <- true.connections %>% filter(view == "para") %>%
      left_join(misty.results_all$importances %>% filter(view == "para.10") , by = c("node1" = "Predictor", "node2" = "Target")) %>%
      filter(!is.na(present), !is.na(Importance))
    
    joined_para_all =joined_para_all %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE)))  
    
    para.roc_all<- joined_para_all  %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    para.pr_all <- joined_para_all %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE))) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    
    # para with pre-filtering CT1
    trim = 0
    targets <- misty.results_CT1$improvements.stats %>% filter(measure == "gain.R2", mean > trim) %>% pull(target)
    
    gain <- misty.results_CT1$improvements.stats %>% filter(measure == "gain.R2")
    
    joined_para_CT1 <- true.connections %>% filter(view == "para") %>%
      left_join(misty.results_CT1$importances %>% filter(view == "para.10") , by = c("node1" = "Predictor", "node2" = "Target")) %>%
      filter(!is.na(present), !is.na(Importance))
    
    joined_para_CT1 =joined_para_CT1 %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE)))  
    
    para.roc_CT1<- joined_para_CT1  %>%
      summarise(auc =  roc.curve(Importance, weights.class0 = present)$auc) %>%
      pull(auc)
    
    para.pr_CT1 <- joined_para_CT1 %>% group_by(sample) %>%
      mutate(Importance = ifelse(node2 %in% targets, Importance, min(Importance,na.rm = TRUE))) %>%
      summarise(auc =  pr.curve(Importance, weights.class0 = present)$auc.integral) %>%
      pull(auc)
    
    
    
    viodata.roc <- data.frame(intra.roc_noCT = intra.roc_all %>% unlist,
                              intra.roc_CT1 = intra.roc_CT1 %>% unlist,
                              para.roc_noCT = para.roc_all %>% unlist,
                              para.roc_CT1 = para.roc_CT1 %>% unlist) %>% 
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.roc, aes(x=Type, y=AUC)) + 
      geom_jitter(aes(col=Type),width = 0.1) + 
      geom_hline(yintercept=0.5, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=0.5, color="#00BFC4", linetype="dashed") +
      theme(axis.text = element_text(family = "Arial")) +
      theme_classic() + ylim(0.0,1)
    
    ggsave(paste0(plot_folder,"AUC_ROC_intra_para.pdf"), width = 4,height = 3.2)
    
    intra.intercept <- sum(true.connections_all %>% filter(view == "intra", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "intra", !is.na(present)) %>% nrow())
    para.intercept <- sum(true.connections_all %>% filter(view == "para", !is.na(present)) %>% pull(present))/
      (true.connections %>% filter(view == "para", !is.na(present)) %>% nrow())
    
    viodata.pr <- data.frame(intra.pr_all = intra.pr_all %>% unlist,
                             para.pr_all = para.pr_all %>% unlist,
                             intra.pr_CT1 = intra.pr_CT1 %>% unlist,
                             para.pr_CT1 = para.pr_CT1 %>% unlist
    ) %>%
      pivot_longer(names_to = "Type", values_to = "AUC", cols = everything())
    
    ggplot(viodata.pr, aes(x=Type, y=AUC)) + 
      geom_jitter(aes(col=Type),width = 0.1) + 
      geom_hline(yintercept=intra.intercept, color="#F8766D", linetype="dashed") + 
      geom_hline(yintercept=para.intercept, color="#00BFC4", linetype="dashed") +
      theme_classic() + ylim(0,1)
    
    ggsave(paste0(plot_folder,"AUC_PR_intra_para.pdf"), width = 4,height = 3.2)
    
    
    # Plot aggregated ROC and PR curves
    
    # for iso-F1 and iso-J lines
    r <- seq(0,1,by=1e-3)
    ks <- c(0.1,0.2,0.5,0.8)
    
    
    
    para.col_all = RColorBrewer::brewer.pal(7,"Paired")[1]
    para.col_CT1 = RColorBrewer::brewer.pal(7,"Paired")[2]
    intra.col_all = RColorBrewer::brewer.pal(7,"Paired")[3]
    intra.col_CT1 = RColorBrewer::brewer.pal(7,"Paired")[4]
    
    
    
    pdf(file = paste0(plot_folder,"ROC_para_intra.pdf"),width = 4,height = 4)
    
    tidy.intra.agg <- misty.results_all$importances.aggregated %>% filter(view == "intra")
    joined.intra.agg <- true.connections_all %>% filter(view == "intra") %>% 
      left_join(tidy.intra.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_intra_all <- roc.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                           curve = TRUE, rand.compute = TRUE)
    plot(roc_intra_all, color=intra.col_all, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE)
    
    # CT1
    
    tidy.intra.agg <- misty.results_CT1$importances.aggregated %>% filter(view == "intra")
    joined.intra.agg <- true.connections_CT1 %>% filter(view == "intra") %>% 
      left_join(tidy.intra.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_intra_CT1 <- roc.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
                               curve = TRUE, rand.compute = TRUE)
    plot(roc_intra_CT1, color=intra.col_CT1, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    
    
    # add para
    tidy.para.agg <- misty.results_all$importances.aggregated %>% filter(view == "para.10") %>%
      mutate(Importance = ifelse(Target %in% targets, Importance, min(Importance,na.rm = TRUE)))
    joined.para.agg <- true.connections_all %>% filter(view == "para") %>% 
      left_join(tidy.para.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_para_all = roc.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                         curve = TRUE, rand.compute = TRUE)
    plot(roc_para_all, color=para.col_all, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    # CT1
    tidy.para.agg <- misty.results_CT1$importances.aggregated %>% filter(view == "para.10") %>%
      mutate(Importance = ifelse(Target %in% targets, Importance, min(Importance,na.rm = TRUE)))
    joined.para.agg <- true.connections_CT1 %>% filter(view == "para") %>% 
      left_join(tidy.para.agg, by = c("node1" = "Predictor", "node2" = "Target")) %>% 
      filter(!is.na(present), !is.na(Importance))
    
    roc_para_CT1 = roc.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
                             curve = TRUE, rand.compute = TRUE)
    plot(roc_para_CT1, color=para.col_CT1, rand.plot = TRUE,main = "ROC curve",auc.main=FALSE, add = TRUE)
    
    text(.4,.05,labels = paste("AUC(para_all)=", round(roc_para_all$auc,digits = 2)),adj = 0,col=para.col_all)
    text(.4,.15,labels = paste("AUC(intra_all)=", round(roc_intra_all$auc,digits = 2)),adj = 0,col=intra.col_all)
    text(.4,.25,labels = paste("AUC(para_CT1)=", round(roc_para_CT1$auc,digits = 2)),adj = 0,col=para.col_CT1)
    text(.4,.35,labels = paste("AUC(intra_CT1)=", round(roc_intra_CT1$auc,digits = 2)),adj = 0,col=intra.col_CT1)
    
    
    ks %>% walk(function(k){
      s <- r + k
      s[s>1] <- NA
      lines(r, s, ylim=c(0,1),xlim=c(0,1), col="gray80")
    })
    dev.off()
    
    
    # 
    # # PR curves
    # para.col = "#EC008C"
    # intra.col = "#00A651"
    # 
    # pdf(file = paste0(plot_folder,"PR_para_intra.pdf"),width = 4,height = 4)
    # pr_intra <- pr.curve(joined.intra.agg %>% pull(Importance), weights.class0 = joined.intra.agg %>% pull(present), 
    #                      curve = TRUE)
    # 
    # plot(pr_intra, color=intra.col, auc.main=FALSE)
    # 
    # 
    # pr_para <- pr.curve(joined.para.agg %>% pull(Importance), weights.class0 = joined.para.agg %>% pull(present), 
    #                     curve = TRUE)
    # plot(pr_para, color=para.col, rand.plot = TRUE, auc.main=FALSE, add = TRUE)
    # 
    # 
    # 
    # 
    # 
    # lines(x = c(0,1), y = c(intra.intercept,intra.intercept), lty = 2, col = intra.col)
    # lines(x = c(0,1), y = c(para.intercept, para.intercept), lty = 2, col = para.col)
    # 
    # 
    # ks %>% walk(function(k){
    #   p <- k*r/(2*r - k)
    #   p[p<=0 | p > 1.5] <- NA
    #   lines(r, p,ylim=c(0,1),xlim=c(0,1), col="gray80")
    # })
    # 
    # text(.4,.8,labels = paste("AUC(para)=", round(pr_para$auc.integral,digits = 3)),adj = 0,col=para.col)
    # text(.4,.7,labels = paste("AUC(intra)=", round(pr_intra$auc.integral,digits = 3)),adj = 0,col=intra.col)
    # 
    # 
    # dev.off()
    # 











