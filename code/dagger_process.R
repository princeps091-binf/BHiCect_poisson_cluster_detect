#DAGGER behaviour
library(data.tree)
library(tidyverse)
library(igraph)
#--------------------------
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
  node_pval<-chr_pval_tbl$emp.pval
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
  alpha_res_l[[as.character(alpha)]]<-tibble(chr=chromo,node=names(which(rejections)),FDR=alpha,emp.pval=node_pval[names(which(rejections))])
  }
  return(do.call(bind_rows,alpha_res_l))
}
res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
#load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/cl_emp_pval_tbl.Rda")
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/peak_fn_shuffle_cl_pval_tbl.Rda")
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/cl_cage_bin_count_hyp_tbl.Rda")
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/cl_cage_bin_count_pval_tbl.Rda")


multi_cage_cl_tbl<-cl_emp_pval_tbl%>%mutate(nbin=as.numeric(unlist(lapply(strsplit(.$cl,split="_"),'[',2))))%>%filter(cage_n>1 & nbin>1)
alpha_seq<-seq(0.001,0.05,length.out=17)
alpha_seq<-c(0.01,0.05)
cage_tss_pval_tbl<-multi_cage_cl_tbl%>%dplyr::rename(emp.pval=cl.emp.pval)
cage_tss_pval_tbl<-cl_hyp_cage_pval_tbl%>%dplyr::rename(emp.pval=hyp.pval)
cage_tss_pval_tbl<-chr_cl_bin_count_pval

chr_res_l<-vector('list',length(unique(cage_tss_pval_tbl$chr)))
names(chr_res_l)<-unique(cage_tss_pval_tbl$chr)
for(chromo in unique(cage_tss_pval_tbl$chr)){
  print(chromo)
  ## Load computed p-values
  chr_top_cl<-unique(unlist(cage_tss_pval_tbl%>%filter(chr==chromo)%>%dplyr::select(cl)))
  chr_pval_tbl<-cage_tss_pval_tbl%>%filter(chr==chromo)
  load(paste0(res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  ## Prune tree to only contain CAGE-containing clusters
  
  cage_set<-unique(c(chr_top_cl,unique(unlist(node_ancestor[chr_top_cl])))) 
  #rebuild corresponding tree
  Prune(chr_bpt, function(x) x$name %in% cage_set)
  chr_res_l[[chromo]]<-  DAGGER_fn(chromo,chr_bpt,chr_pval_tbl,alpha_seq)
  
}
chr_dagger_tbl<-do.call(bind_rows,chr_res_l)
chr_dagger_tbl<-chr_dagger_tbl%>%left_join(.,cage_tss_pval_tbl%>%dplyr::select(chr,cl,res,emp.pval)%>%dplyr::rename(node=cl))
save(chr_dagger_tbl,file='~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/dagger_naive_fdr_0105_cage_bin_boot_tbl.Rda')
load('~/Documents/multires_bhicect/weeklies/weekly29/chr_dagger_tbl.Rda')
chr_dagger_tbl%>%group_by(res)%>%summarise(n=n())%>%ggplot(.,aes(res,n))+geom_bar(stat='identity')
gg_chr<-chr_dagger_tbl%>%group_by(chr)%>%distinct(node)%>%summarise(n=n())%>%ggplot(.,aes(chr,n))+geom_bar(stat='identity')
gg_chr<-gg_chr + theme_classic()
gg_chr
ggsave("~/Documents/multires_bhicect/weeklies/weekly29/img/dagger_chr.png",gg_chr)
gg_pval<-chr_dagger_tbl%>%ggplot(.,aes(-log10(emp.pval)))+geom_histogram(bins=60)+facet_wrap(FDR~.,scales="free")
gg_pval<-gg_pval + theme_classic()
gg_pval
ggsave("~/Documents/multires_bhicect/weeklies/weekly29/img/dagger_pval.png",gg_pval)
#-------------------------------------------------------------------------------------------------------
library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
dagger_set_tbl<-chr_dagger_tbl%>%distinct(chr,nodes,emp.pval)%>%dplyr::rename(cl=nodes)%>%left_join(.,cage_tss_pval_tbl)
#load the CAGE annotation
load("~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_coord_tbl.Rda")
cage_hmec_a<-cage_hmec_a%>%filter(!(is.na(start)))
#Build the corresponding GRanges object to compute overlap
cage_hmec_Grange<-   GRanges(seqnames=cage_hmec_a$chr,
                             ranges = IRanges(start=cage_hmec_a$start,
                                              end=cage_hmec_a$end,
                                              names=paste(cage_hmec_a$chr,1:nrow(cage_hmec_a),sep='_')
                             ))
mcols(cage_hmec_Grange)<-tibble(namess=cage_hmec_a$namess,gene=cage_hmec_a$ENSG,I=cage_hmec_a$m)
top_dagger_Grange<-c()
for(chromo in unique(dagger_set_tbl$chr)){
  print(chromo)
  
  load(paste0("~/Documents/multires_bhicect/data/HMEC/spec_res/",chromo,"_spec_res.Rda"))
  #chr_cl_tbl<-cage_union_pval_tbl%>%filter(chr==chromo & emp.pval<=0.0002)
  chr_cl_tbl<-dagger_set_tbl%>%filter(chr==chromo)
  chr_cl_tbl<-chr_cl_tbl%>%mutate(bins=chr_spec_res$cl_member[chr_cl_tbl$cl])
  tmp_cl_tbl<-chr_cl_tbl%>%unnest(bins)%>%mutate(start=as.numeric(bins))%>%mutate(end=start + res_num[res]-1)%>%dplyr::select(-bins)
  
  cl_Grange <- GRanges(seqnames=tmp_cl_tbl$chr,
                       ranges = IRanges(start=tmp_cl_tbl$start,
                                        end=tmp_cl_tbl$end,
                                        names=paste(tmp_cl_tbl$chr,1:nrow(tmp_cl_tbl),sep='_')
                       ))
  cl_Grange<-IRanges::reduce(cl_Grange)
  top_dagger_Grange<-append(top_dagger_Grange,cl_Grange)
  

}
#--------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(formattable)

dagger_gene_content<-unique(unlist(cage_hmec_Grange@elementMetadata$gene[unique(subjectHits(findOverlaps(top_dagger_Grange,cage_hmec_Grange)))]))

cage_active_gene<-unique(unlist(cage_hmec_a$ENSG))
ego_n <- enrichGO(gene        = dagger_gene_content,
                  universe      = cage_active_gene,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.001,
                  qvalueCutoff  = 0.001,
                  readable      = TRUE)
formattable(ego_n@result%>%dplyr::select(Description,p.adjust)%>%arrange(p.adjust)%>%dplyr::slice_head(n=13))

library(biomaRt)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

ensembl <- useMart("ensembl")
#listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
#listAttributes(ensembl)


hm_gene_set <- read_delim("~/Documents/multires_bhicect/data/epi_data/Gene_annotation/c2.all.v7.3.entrez.gmt",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)

hallmark_set<-lapply(1:nrow(hm_gene_set),function(x){
  tmp<-hm_gene_set[x,-c(1,2)]
  return(tmp[!(is.na(tmp))])})
names(hallmark_set)<-hm_gene_set$X1

genes_hm <- getBM(
  filters="entrezgene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=unique(unlist(hallmark_set)),
  mart=ensembl)


genes_hmec <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=cage_active_gene,
  mart=ensembl)

hg19_ENST_to_ENSG <- read_delim("~/Documents/multires_bhicect/data/epi_data/hg19_ENST_to_ENSG.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
hg19_gene<-unique(hg19_ENST_to_ENSG$name2)

genes_hg19 <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=hg19_gene,
  mart=ensembl)

tot_gene_conv<-tibble(genes_hmec%>%full_join(.,genes_hm)%>%full_join(.,genes_hg19))

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
})
clusterExport(cl,c('tot_gene_conv','dagger_gene_content','cage_active_gene'))
go_pval<-parLapply(cl,hallmark_set,function(tmp_set){
  hitInSample<-nrow(tot_gene_conv%>%filter(ensembl_gene_id %in% dagger_gene_content)%>%filter(entrezgene_id %in% tmp_set)%>%distinct(entrezgene_id))
  sampleSize<-nrow(tot_gene_conv%>%filter(ensembl_gene_id %in% dagger_gene_content)%>%distinct(entrezgene_id))
  hitInPop<-nrow(tot_gene_conv%>%filter(ensembl_gene_id %in% cage_active_gene)%>%filter(entrezgene_id %in% tmp_set)%>%distinct(entrezgene_id))
  failInPop<-nrow(tot_gene_conv%>%filter(ensembl_gene_id %in% cage_active_gene)%>%distinct(entrezgene_id)) - hitInPop
  return(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
})
stopCluster(cl)
rm(cl)
top_set<-head(sort(p.adjust(unlist(go_pval),method='fdr')),13)
path_tbl<-tibble(Gene.Set=names(top_set),p.adjust=top_set)
formattable(path_tbl%>%arrange(p.adjust))
#---------------------
"1Mb_1_0_87000000_87000000"