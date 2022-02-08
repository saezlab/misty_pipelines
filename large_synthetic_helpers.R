


get_connections <- function(current_cell_type, true.interactions, all_markers) {
    
    # Intra
    true.intra.connections <- get_intra_connections(current_cell_type,true.interactions)
        
    # Intercellular signaling. 
    all_inter.connections <- get_inter_connections(current_cell_type,true.interactions,all_markers)
    
    
    true.connections = bind_rows(true.intra.connections,all_inter.connections)
}

# Intra
# ground truth for intracellular signaling: 
# - interactions that represents the ground truth are markerd as present = TRUE. 
# - to get a proper AUC curves we have to enumerate all the possible cases
#
# case of "all" cell types: 
# - we merge the possible cases (across all cell types)
# we remove the interactions that are with Ligands because ligand measurements 
# are removed from the data
#
# case of cellType X:
#- we filter the interactions that exist in cell type X. 
#- remove interactions with ligands because ligand measurements were removed from data. 
get_intra_connections <- function(current_cell_type,true.interactions){
    
    if(current_cell_type == "all"){
        true.intra.connections <- true.interactions %>% rename(node1 = "from",node2="target") %>%
            mutate(present = Vmax>0,
                   view = "intra") %>%
            mutate(cell_type = "all") %>%
            # remove ligands: they are also removed from misty inputs
            filter(!grepl("L.",node1)) %>%
            filter(!grepl("L.",node2)) %>%
            select(node1,node2, present,view,cell_type) %>% 
            # remove duplicates that appeared in more cell types
            unique() %>%
            mutate(depth = 1)
        
    }else{
        true.intra.connections <- true.interactions %>% rename(node1 = "from",node2="target") %>%
            mutate(present = Vmax>0,
                   view = "intra") %>%
            filter(!grepl("L.",node1)) %>%
            filter(!grepl("L.",node2)) %>%
            filter(cell_type==current_cell_type) %>%
            select(node1,node2, present,view,cell_type) %>%
            mutate(depth = 1)
    }
    return(true.intra.connections)
}


get_inter_connections <- function(current_cell_type,true.interactions,all_markers){
    
    
    # Intercellular signaling. 
    # all possible interactions: all pairs of markers
    
    # ground truth: 
    # def: there is a cell-cell communication between two cell types (A and B) through ligand_X if CT_A  expresses
    # ligand_X and CT_B has receptor_X which detects that. 
    
    intra <- true.interactions %>% mutate(present = Vmax > 0) %>% 
        filter(present) %>%
        select(from,target,cell_type)
    
    # Ligand production side -------
    # X1 -> L
    ld1 <- intra %>% 
        filter(grepl("^L",target)) %>%
        mutate(depth = 1)
    
    # X2 -> X1 -> L   =>    X2 -> L1
    ld2 <- ld1 %>% 
        inner_join(intra,by=c("from"="target", "cell_type"="cell_type")) %>%
        # we drop the intermediate node
        select(-from) %>%
        rename(from = "from.y") %>%
        select(from, target, cell_type) %>%
        mutate(depth = 2)
    
    
    # X3 -> X2 -> X1 -> L   =>    X3 -> L1
    ld3 <- ld2 %>% 
        inner_join(intra,by=c("from"="target", "cell_type"="cell_type")) %>%
        # we drop the intermediate node
        select(-from) %>%
        rename(from = "from.y") %>%
        select(from, target, cell_type) %>%
        mutate(depth = 3)
    
    
    # Receptor production side -------
    # D0: L - > R
    rd0 <- intra %>% 
        filter(grepl("^R",target)) %>%
        mutate(depth = 0)
    
    
    # R -> X1
    rd1 <- rd0 %>% 
        inner_join(intra,by=c("target"="from", "cell_type"="cell_type")) %>%
        # we drop the intermediate node
        select(-target) %>%
        rename(target = "target.y") %>%
        select(from, target, cell_type) %>%
        mutate(depth = 1)
    
    # R -> X1 -> X2   =>    R -> X2
    rd2 <- rd1 %>% 
        inner_join(intra,by=c("target"="from", "cell_type"="cell_type")) %>%
        # we drop the intermediate node
        select(-target) %>%
        rename(target = "target.y") %>%
        select(from, target, cell_type) %>%
        mutate(depth = 2)
    
    
    # R -> X1 -> X2 ->X3   =>    R -> X3
    rd3 <- rd2 %>% 
        inner_join(intra,by=c("target"="from", "cell_type"="cell_type")) %>%
        # we drop the intermediate node
        select(-target) %>%
        rename(target = "target.y") %>%
        select(from, target, cell_type) %>%
        mutate(depth = 3)
    
    
    ligand_production = bind_rows(ld1,ld2,ld3)
    ligand_uptake = bind_rows(rd0,rd1,rd2,rd3)
    
    
    # combine production and uptake such that the target of production (ligandX) should
    # match the activator in the uptake (ligandX). 
    inter.connections_base <- full_join(ligand_production, ligand_uptake, by=c("target"="from"),suffix = c("_source", "_target")) %>%
        # the target (ligand) now is called via (the mediator) and 
        # target of the ligand in the target cell type (target_target) is simplified to target. 
        rename(via = "target",
               target = "target_target") %>%
        mutate(depth = depth_target + depth_source-1) %>%
        select(from,  target, cell_type_source,  cell_type_target, depth)
    
    # if we discard cell-type information, then we should remove this information and 
    # take the unique interactions. 
    # if we consider cell types, we consider only the type of the target cell, so 
    # we drop the cell type info on the source cell and filter for the current cell type. 
    if(current_cell_type == "all"){
        true.inter.connections <- inter.connections_base %>%
            select(-cell_type_source, -cell_type_target) %>%
            unique() %>%
            mutate(cell_type = "all") 
    }else{
        true.inter.connections <- inter.connections_base %>%
            select(-cell_type_source) %>%
            filter(cell_type_target == current_cell_type) %>% unique() %>%
            rename(cell_type = cell_type_target)
    }
    
    # so far we have the true interactions, but we also need to enumerate all possible interactions
    # this interactome is all possible pairs of species. 
    
    all_possible_intercellular_interactions <- expand.grid(all_markers,all_markers) %>%
        rename(from = Var1,target = Var2) %>% as_tibble() %>%
        mutate(cell_type = current_cell_type)
    
    # then we merge the two TRUE and all combinations:
    all_inter.connections <- true.inter.connections %>% 
         unique() %>%
        mutate(present = TRUE) %>% 
        full_join(all_possible_intercellular_interactions, by=c("from","target","cell_type")) %>%
        mutate(present=ifelse(is.na(present),FALSE,present)) %>%
        mutate(present = ifelse(from==target,NA,present)) %>%
        mutate(view = "para") %>%
        rename(node1 = from,
               node2 = target)
    return(all_inter.connections)   
}

