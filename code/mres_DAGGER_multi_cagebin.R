#DAGGER stratified by resolution
library(data.tree)
library(tidyverse)
library(igraph)
#--------------------------
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
rej_fn<-function(nodes,lvl,num_rejected,alpha){
  #pick node effective node and leaf number
  ms_d<-ms[nodes]
  ls_d<-ls[nodes]
  p_vals_d<- node_pval[nodes]
  
  
  ### P-value threshold function
  # r is the considered rank
  crit_func <- function(r,alpha){
    alpha * ls_d * (ms_d + r + num_rejected[as.character(lvl-1)] - 1) / l / ms_d
  }
  
  r <- length(p_vals_d)
  # Determine the appropriate threshold value
  while (sum(p_vals_d <= crit_func(r,alpha)) < r){
    r <- r-1 
  }  
  R <- r 
  tmp_rejected<-nodes[which(p_vals_d <= crit_func(R,alpha))]
  
  return(tmp_rejected)
}

detect_inter_cage_cl_fn<-function(feature_coord_tbl,feature_pval_tbl,res_num){
  
  fn_env<-environment()
  
  feature_Grange<-   GRanges(seqnames=feature_coord_tbl$chr,
                             ranges = IRanges(start=feature_coord_tbl$start,
                                              end=feature_coord_tbl$end
                             ))
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("feature_pval_tbl","feature_Grange","res_num"),envir = fn_env)
  feature_bin_n_l<-parLapply(cl,1:nrow(feature_pval_tbl),function(x){
    cl_Grange<-   GRanges(seqnames=feature_pval_tbl$chr[x],
                          ranges = IRanges(start=as.numeric(feature_pval_tbl$bins[[x]]),
                                           end=as.numeric(feature_pval_tbl$bins[[x]]) + res_num[feature_pval_tbl$res[x]]-1
                          ))
    return(length(unique(queryHits(findOverlaps(cl_Grange,feature_Grange)))))
    
  })
  stopCluster(cl)
  rm(cl)
  feature_pval_tbl<-feature_pval_tbl %>% mutate(feature.bin=unlist(feature_bin_n_l))
  return(feature_pval_tbl)
  
}

#--------------------------
#load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/cl_emp_pval_tbl.Rda")
load("./data/pval_tbl/CAGE_tss_HMEC_pval_tbl.Rda")
cage_tss_pval_tbl<-cl_chr_emp_pval_tbl
load("./data/CAGE_tss_coord_HMEC_tbl.Rda")
cage_coord_tbl<-cage_tss_HMEC_tbl%>%filter(!(is.na(start)))
res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"


chr_res_l<-vector('list',length(unique(cage_tss_pval_tbl$chr)))
names(chr_res_l)<-unique(cage_tss_pval_tbl$chr)
alpha_seq<-0.01

for(chromo in unique(cage_tss_pval_tbl$chr)){
  print(chromo)
  # Build the BHiCect tree
  load(paste0(res_file,chromo,"_spec_res.Rda"))

  chr_pval_tbl<-cage_tss_pval_tbl %>% filter(chr==chromo)
  # select cluster with at least two CAGE-containing bins (trx aggregation/looping)
  chr_pval_tbl<-detect_inter_cage_cl_fn(cage_coord_tbl,chr_pval_tbl,res_num) %>% 
    filter(feature.bin>1)
  
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  #collect leaves and ancestors
  tmp_leaves<-chr_bpt$Get('name',filterFun=isLeaf)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  #build the cage-containing sub-tree
  cage_node<-unlist(chr_pval_tbl%>%dplyr::select(cl))
  cage_set<-unique(c(cage_node,unique(unlist(node_ancestor[cage_node])))) 
  Prune(chr_bpt, function(x) x$name %in% cage_set)
  p_node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  
  #rebuild corresponding tree according to DAGGER
  node_dagger_children<-lapply(p_node_ancestor,'[',2)
  #eleminate Root node to "create" DAGGER leaves
  node_dagger_children<-node_dagger_children[-1]
  
  node_dagger_children<-lapply(node_dagger_children,function(x){
    if(x == "Root"){return(NULL)} else{
      return(x)
    }
  })
  #Create DAGGER-parent mapping (immediate BPT children of each node)
  node_dagger_parent<-lapply(cage_set,function(x){
    tmp<-names(which(unlist(node_dagger_children) == x))
    return(unlist(lapply(strsplit(tmp,split="\\."),'[',1)))
  })
  names(node_dagger_parent)<-cage_set
  
  g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')
  
  tmp_res_l<-vector('list',length(unique(chr_pval_tbl$res)))
  names(tmp_res_l)<-unique(chr_pval_tbl$res)
  for(tmp_res in unique(chr_pval_tbl$res)){
    print(tmp_res)
    res_cage_node<-unlist(chr_pval_tbl%>%filter(res==tmp_res)%>%dplyr::select(cl))
    res_cage_set<-unique(c(res_cage_node,grep(paste0(tmp_res),unique(unlist(node_ancestor[res_cage_node])),value=T))) 
    if(length(res_cage_set)<2){
      tmp_pval<-unlist(chr_pval_tbl%>%filter(cl%in%res_cage_set)%>%dplyr::select(emp.pval))
      tmp_res_l[[tmp_res]]<-tibble(chr=chromo,res=tmp_res,node=res_cage_set,FDR=NA,emp.pval=tmp_pval)
      next
    }
    # Produce CAGE-tree for each resolution
    tmp_g<-induced_subgraph(g_bpt,res_cage_set)
    tmp_g_comp<-components(tmp_g)
    ## detect the DAGGER leaves
    ### For each component detect top parent node and label it as leaf
    #### DAGGER leaves are nodes who don't have any BPT-ancestors among the node at the considered resolution! 
    dagger_leaf<-names(which(unlist(lapply(node_dagger_children[res_cage_set],function(x)sum(grepl(tmp_res,x))))<1))
    
    dagger_roots<-res_cage_set[!(res_cage_set %in% unique(unlist(lapply(p_node_ancestor[res_cage_set],'[',-1))))]
    #assign levels/depth to this node set
    ## For each component, depth is the graph distance with the corresponding DAGGER-roots (BPT-leaves)
    comp_set<-which(tmp_g_comp$csize > 1)
    comp_tbl<-do.call(bind_rows,lapply(comp_set,function(x){
      tmp_node<-names(which(tmp_g_comp$membership==x))
      #not working
      tmp_d<-distances(tmp_g,tmp_node[which(tmp_node %in% dagger_roots)],tmp_node,mode = 'in')
      tmp_lvl<-apply(tmp_d,2,function(x)max(x[!(is.infinite(x))])+1)
      tmp_lvl[names(which(apply(tmp_d,2,function(a)any(a==0))))]<-1
      return(tibble(lvl=tmp_lvl,node=names(tmp_lvl),comp=x))
    }))
    gen_lvl<-max(comp_tbl$lvl)
    if(is.infinite(gen_lvl)){gen_lvl<-1}
    isl_set<-which(tmp_g_comp$csize == 1)
    isl_tbl<-do.call(bind_rows,lapply(isl_set,function(x){
      tmp_node<-names(which(tmp_g_comp$membership==x))
      #not working
      return(tibble(lvl=gen_lvl,node=tmp_node,comp=x))
    }))
    node_m_lvl_tbl<-comp_tbl%>%bind_rows(.,isl_tbl)%>%dplyr::rename(m_lvl=lvl)
    #-------------------------------------
    l<-length(dagger_leaf)
    #Build the node p-value mapping vector
    node_pval<-unlist(chr_pval_tbl%>%filter(cl %in% res_cage_set)%>%dplyr::select(emp.pval))
    names(node_pval)<-unlist(chr_pval_tbl%>%filter(cl %in% res_cage_set)%>%dplyr::select(cl))
    
    #Compute recursively the effective number of leaves and nodes
    ### Computation for effective leaf/node number is done from DAGGER-leaves to DAGGER-roots
    ms<-rep(0,nrow(node_m_lvl_tbl))
    names(ms)<-node_m_lvl_tbl$node
    ls<-rep(0,nrow(node_m_lvl_tbl))
    names(ls)<-node_m_lvl_tbl$node
    #assign 1 to DAGGER-leaf clusters
    ls[dagger_leaf]<-1
    ms[dagger_leaf]<-1
    #Loop through levels in decreasing DAGGER depth order (starting with DAGGER-leaves)
    for(lvl in rev(sort(unique(node_m_lvl_tbl$m_lvl)))){
      #print(lvl)
      tmp_node<-unlist(node_m_lvl_tbl%>%filter(m_lvl==lvl)%>%dplyr::select(node))
      tmp_leaves_idx<-which(tmp_node %in% dagger_leaf)
      if(length(tmp_leaves_idx)>0){tmp_node<-tmp_node[-tmp_leaves_idx]}
      if(length(tmp_node)<1){next}
      for(p in tmp_node){
        
        ms[p]<-1+sum(ms[node_dagger_children[[p]]]/unlist(lapply(node_dagger_children[[p]],function(x)length(node_dagger_parent[[x]]))))
        
        
        ls[p]<-sum(ls[node_dagger_children[[p]]]/unlist(lapply(node_dagger_children[[p]],function(x)length(node_dagger_parent[[x]]))))
      }
    }
    #Filter out clusters with too few bins or CAGE peaks
    node_m_lvl_tbl<-node_m_lvl_tbl%>%mutate(nbin=as.numeric(unlist(lapply(strsplit(.$node,split="_"),'[',2))))%>%filter(nbin>1)
    alpha_res_l<-vector('list',length(alpha_seq))
    names(alpha_res_l)<-as.character(alpha_seq)
    for ( alpha in alpha_seq){
      
      #----------------------------------------
      print("Node rejection")
      #vector recording the number of rejections at every depth
      num_rejected = rep(0, 1 + max(node_m_lvl_tbl$m_lvl)) 
      names(num_rejected)<-as.character(seq(0, max(node_m_lvl_tbl$m_lvl)))
      #vector recording the actual nodes rejecting the null hypothesis
      rejections<-rep(F,nrow(node_m_lvl_tbl))
      names(rejections)<-node_m_lvl_tbl$node
      #loop through increasing depth levels (starting from DAGGER-roots)
      for (lvl in seq(1, max(node_m_lvl_tbl$m_lvl))){
        #print(lvl)
        nodes_depth_d <- unlist(node_m_lvl_tbl%>%filter(m_lvl==lvl)%>%dplyr::select(node)) 
        
        # Delete the nodes one of whose parents has not been rejected.
        if ( lvl > 1){
          
          nodes_depth_d<-nodes_depth_d[unlist(lapply(nodes_depth_d,function(x){
            all(node_dagger_parent[[x]] %in% names(which(rejections)))
          }))]
          if(any(is.na(node_pval[nodes_depth_d]))){
            #further filter out the nodes for which we don't have p-values
            nodes_depth_d<-names(which(!(is.na(node_pval[nodes_depth_d]))))
            
          }
          
        }
        # Performs the rejection step at depth d.  
        rejected_nodes_depth_d <- rej_fn(nodes_depth_d, lvl, num_rejected,alpha)
        rejections[rejected_nodes_depth_d] <- T
        num_rejected[as.character(lvl)] <- num_rejected[as.character(lvl-1)] + length(rejected_nodes_depth_d)
        
      }
      alpha_res_l[[as.character(alpha)]]<-tibble(chr=chromo,res=tmp_res,node=names(which(rejections)),FDR=alpha,emp.pval=node_pval[names(which(rejections))])
      
    }
    tmp_res_l[[tmp_res]]<-do.call(bind_rows,alpha_res_l)
    
  }
  
  chr_res_l[[chromo]]<-do.call(bind_rows,tmp_res_l)
}

dagger_mres_tbl<-do.call(bind_rows,chr_res_l)
dagger_mres_tbl%>%ggplot(.,aes(emp.pval))+geom_histogram()+facet_wrap(FDR~.,scales="free")
dagger_mres_tbl<-dagger_mres_tbl%>%filter(!(is.na(FDR)))
save(dagger_mres_tbl,file="./data/DAGGER_tbl/HMEC_tss_dagger_tbl.Rda")

dagger_mres_tbl2<-dagger_mres_tbl
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/dagger_mres_fdr_01_multi_cagebin_tbl.Rda")
#Some level of discrepency -> need to evaluate consistency in gene-set enrichment
dagger_mres_tbl2 %>% dplyr::select(chr,node) %>% left_join(.,dagger_mres_tbl%>% dplyr::select(chr,node,FDR)) %>% filter(is.na(FDR))