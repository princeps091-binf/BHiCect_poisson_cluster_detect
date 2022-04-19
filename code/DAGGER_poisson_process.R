library(GenomicRanges)
library(data.tree)
library(tidyverse)
library(igraph)
#--------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
obj_in_fn<-function(file){
  obj_out<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(obj_out)
}

detect_inter_cage_cl_fn<-function(feature_GRange,feature_pval_tbl,chr_spec_res,res_num){
  
  fn_env<-environment()
  
  feature_pval_tbl<-feature_pval_tbl %>% 
    mutate(bins=chr_spec_res$cl_member[feature_pval_tbl$cl])
  

  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("feature_pval_tbl","feature_GRange","res_num"),envir = fn_env)
  feature_bin_n_l<-parLapply(cl,1:nrow(feature_pval_tbl),function(x){
    cl_Grange<-   GRanges(seqnames=feature_pval_tbl$chr[x],
                          ranges = IRanges(start=as.numeric(feature_pval_tbl$bins[[x]]),
                                           end=as.numeric(feature_pval_tbl$bins[[x]]) + res_num[feature_pval_tbl$res[x]]-1
                          ))
    return(length(unique(queryHits(findOverlaps(cl_Grange,feature_GRange)))))
    
  })
  stopCluster(cl)
  rm(cl)
  feature_pval_tbl<-feature_pval_tbl %>% mutate(feature.bin=unlist(feature_bin_n_l))
  return(feature_pval_tbl)
  
}

DAGGER_fn<-function(chromo,chr_bpt,chr_pval_tbl,alpha){
  print("Tree building")
  cage_set<-chr_bpt$Get('name')
  p_node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  #p_node_ancestor<-lapply(p_node_ancestor,'[',-1)
  cage_bpt_leaf<-chr_bpt$Get('name',filterFun = isLeaf)
  node_bpt_lvl<-chr_bpt$Get('level')
  #----------------------------------------
  print("DAGGER genealogy building")
  
  # Create DAGGER-children mapping ~ immediate BPT-parent
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
  # Assign deepest possible level for each node wrt connected bpt-leaf
  bpt_lvl_track<-do.call(bind_rows,lapply(cage_bpt_leaf,function(x){
    tmp_vec<-p_node_ancestor[[x]]
    tmp<-rev(node_bpt_lvl[tmp_vec])
    names(tmp)<-tmp_vec
    
    return(tibble(node=tmp_vec,lvl=tmp,leaf=x)%>%dplyr::slice(n=-n()))}))
  node_m_lvl_tbl<-bpt_lvl_track%>%bind_rows(.,tibble(node=cage_bpt_leaf,lvl=1,leaf=cage_bpt_leaf))%>%group_by(node)%>%summarise(m_lvl=max(lvl))
  #----------------------------------------
  dagger_leaf<-names(which(unlist(lapply(node_dagger_children,length))<1))
  
  l<-length(dagger_leaf)
  
  #----------------------------------------
  print("Build node vectors")
  
  #Build the node p-value mapping vector
  node_pval<-chr_pval_tbl$pois.pval
  names(node_pval)<-chr_pval_tbl$cl
  
  ## Compute the effective number of nodes and leaves for each node
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
  #eliminate nodes with less than 2 constitutive bins
  node_m_lvl_tbl<-node_m_lvl_tbl%>%mutate(nbin=as.numeric(unlist(lapply(strsplit(.$node,split="_"),'[',2))))%>%filter(nbin>1)
  
  #----------------------------------------
  rej_fn<-function(nodes,node_pval,lvl,num_rejected,alpha){
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
  #----------------------------------------
  alpha_res_l<-vector('list',length(alpha_seq))
  names(alpha_res_l)<-as.character(alpha_seq)
  for ( alpha in alpha_seq){
    
    print(paste("Node rejection:",alpha))
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
      rejected_nodes_depth_d <- rej_fn(nodes_depth_d,node_pval, lvl, num_rejected,alpha)
      rejections[rejected_nodes_depth_d] <- T
      num_rejected[as.character(lvl)] <- num_rejected[as.character(lvl-1)] + length(rejected_nodes_depth_d)
      
    }
    alpha_res_l[[as.character(alpha)]]<-tibble(chr=chromo,node=names(which(rejections)),FDR=alpha,pois.pval=node_pval[names(which(rejections))])
  }
  return(do.call(bind_rows,alpha_res_l))
}

#--------------------------
res_folder<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
pval_tbl_file<-"./data/pval_tbl/HMEC_pois_pval_tbl.Rda"
feature_GRange_file<-"./data/GRanges/CAGE_union_HMEC_Grange.Rda"
out_file<-"./data/pval_tbl/DAGGER/HMEC_poisson_DAGGER_05.Rda"
#--------------------------
pval_tbl<-obj_in_fn(pval_tbl_file)
feature_GRange<-obj_in_fn(feature_GRange_file)

chr_set<-unique(pval_tbl$chr)
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set

for(chromo in chr_set){
  message(chromo)
  # Build the BHiCect tree
  chr_spec_res<-obj_in_fn(paste0(res_folder,chromo,"_spec_res.Rda"))

  chr_pval_tbl<-pval_tbl %>% filter(chr==chromo)
  
  # select cluster with at least two CAGE-containing bins (trx aggregation/looping)
  chr_pval_tbl<-detect_inter_cage_cl_fn(feature_GRange,chr_pval_tbl,chr_spec_res,res_num) %>% 
    filter(feature.bin>1)
  
  chr_top_cl<-unique(unlist(chr_pval_tbl%>%dplyr::select(cl)))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  ## Prune tree to only contain CAGE-containing clusters
  
  cage_set<-unique(c(chr_top_cl,unique(unlist(node_ancestor[chr_top_cl])))) 
  #rebuild corresponding tree
  Prune(chr_bpt, function(x) x$name %in% cage_set)
  alpha_seq<-0.05
  chr_res_l[[chromo]]<-  DAGGER_fn(chromo,chr_bpt,chr_pval_tbl,alpha_seq) %>% 
    left_join(.,chr_pval_tbl %>% 
                dplyr::select(chr,res,cl,feature_n,pois.pval) %>% 
                dplyr::rename(node=cl))
  
}
chr_dagger_tbl<-do.call(bind_rows,chr_res_l)
save(chr_dagger_tbl,file=out_file)
