library(tidyverse)
library(vroom)
library(parallel)
#---------------------------------------------------------------------------
## Utility functions
cage_tbl_subset<-function(cage_tbl,ID_col,col_set){
  return(cage_tbl%>%dplyr::select(ID_col)%>%bind_cols(cage_tbl%>%dplyr::select(contains(col_set))))
}
cage_sub_fn<-function(cage_tbl){
  print("compute m")
  cl<-makeCluster(5)
  tmp_m<-parallel::parApply(cl,X = as.matrix(cage_tbl[,-1]),MARGIN = 1,function(x){
    mean(x)
  })
  stopCluster(cl)
  rm(cl)
  cage_tbl<-cage_tbl%>%mutate(m=tmp_m)%>%filter(m>0)
  return(cage_tbl)
}
cage_tbl_coord_build_fn<-function(cage_tbl,ID_col){
  
  cage_coord<-cage_tbl %>% dplyr::select(ID_col) %>% unlist
  cage_start<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',1)))
  cage_end<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',3)))
  cage_chr<-unlist(lapply(strsplit(cage_coord,split = ':'),'[',1))
  return(cage_tbl%>%mutate(chr=cage_chr,start=cage_start,end=cage_end))
  
}

#---------------------------------------------------------------------------

# CAGE-peak annotation
base::load("~/Documents/multires_bhicect/data/epi_data/CAGE/CAGE_tss_ann_smpl.Rda")
# CAGE-peak for TSS with tpm normalisation
cage<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#')
colnames(cage)[1]<-"Id"
## raw count TSS data
cage_count<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_counts_ann.osc.txt",comment = '#')

# Sample set corresponding to the cell-line of interest
## grep("K562",colnames(cage),value = T)

MDA<-"CNhs10736"
K562<-c("CNhs12334", "CNhs12335", 'CNhs12336', 'CNhs11250', 'CNhs12458', 'CNhs12684','CNhs12786')
MCF7<-c("CNhs11943","CNhs12564","CNhs12475","CNhs12703")

H1<-c("CNhs14067","CNhs14068","CNhs13964")
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')

cage_tss_H1_tbl<-cage %>% cage_tbl_subset(.,"Id",H1) %>%
  cage_sub_fn()
cage_tss_HMEC_tbl<-cage %>% cage_tbl_subset(.,"Id",HMEC) %>%
  cage_sub_fn()
cage_tss_GM12878_tbl<-cage %>% cage_tbl_subset(.,"Id",GM12878) %>%
  cage_sub_fn()

cage_tss_H1_tbl<-cage_tss_H1_tbl %>% cage_tbl_coord_build_fn(.,"Id")
cage_tss_HMEC_tbl<-cage_tss_HMEC_tbl %>% cage_tbl_coord_build_fn(.,"Id")
cage_tss_GM12878_tbl<-cage_tss_GM12878_tbl %>% cage_tbl_coord_build_fn(.,"Id") 

cage_tss_H1_tbl<-cage_tss_H1_tbl %>% left_join(.,tss_ann,by=c("Id"="namess"))
cage_tss_HMEC_tbl<-cage_tss_HMEC_tbl %>% left_join(.,tss_ann,by=c("Id"="namess"))
cage_tss_GM12878_tbl<-cage_tss_GM12878_tbl %>% left_join(.,tss_ann,by=c("Id"="namess"))

save(cage_tss_H1_tbl,file="./data/CAGE_tss_coord_H1_tbl.Rda")
save(cage_tss_HMEC_tbl,file='./data/CAGE_tss_coord_HMEC_tbl.Rda')
save(cage_tss_GM12878_tbl,file="./data/CAGE_tss_coord_GM12878_tbl.Rda")
