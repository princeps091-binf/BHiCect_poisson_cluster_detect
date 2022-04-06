#Get cell-line specific cage-data and active cage-clusters
library(tidyverse)
library(vroom)
library(parallel)

# CAGE-peak for enhancers
cage_enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",comment = '#')

#---------------------------------------------------------------------------
## Utility functions
cage_tbl_subset<-function(cage_tbl,ID_col,col_set){
  return(cage_tbl%>%dplyr::select(ID_col)%>%bind_cols(cage_tbl%>%dplyr::select(contains(col_set))))
}
cage_sub_fn<-function(cage,col_set){
  cage_set<-cage%>%dplyr::select(`00Annotation`)%>%bind_cols(cage%>%dplyr::select(contains(col_set)))
  print("compute m")
  cl<-makeCluster(5)
  tmp_m<-parallel::parApply(cl,X = as.matrix(cage_set[,-1]),MARGIN = 1,function(x){
    mean(x)
  })
  stopCluster(cl)
  rm(cl)
  cage_set<-cage_set%>%mutate(m=tmp_m)%>%filter(m>0)
  return(cage_set)
}
cage_tbl_coord_build_fn<-function(cage_tbl,ID_col){
  
  cage_coord<-cage_tbl %>% dplyr::select(ID_col) %>% unlist
  cage_start<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split='-'),'[',1)))
  cage_end<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split='-'),'[',2)))
  cage_chr<-unlist(lapply(strsplit(cage_coord,split = ':'),'[',1))
  return(cage_tbl%>%mutate(chr=cage_chr,start=cage_start,end=cage_end))
  
}

tss_enh_tbl_build_fn<-function(enh_file,tss_file){
  enh_tbl<-get(base::load(enh_file))
  tmp_obj<-names(mget(base::load(enh_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  tss_tbl<-get(base::load(tss_file))
  tmp_obj<-names(mget(base::load(tss_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tss_tbl %>% filter(!is.na(start)) %>% dplyr::select(Id,chr,start,end) %>% 
           bind_rows(.,enh_tbl %>% dplyr::select(Id,chr,start,end)))
  
}

#---------------------------------------------------------------------------

# Identify the set of samples coinciding with the desired cell-type
## grep("K562",colnames(cage),value = T)
###as keyword MDA ~ oestrogen negative

H1<-c("CNhs14067","CNhs14068","CNhs13964")
MDA<-"CNhs10736"
K562<-c("CNhs12334", "CNhs12335", 'CNhs12336', 'CNhs11250', 'CNhs12458', 'CNhs12684','CNhs12786')
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')
MCF7<-c("CNhs11943","CNhs12564","CNhs12475","CNhs12703")


cage_enh_H1_tbl<-cage_enh_tbl %>% cage_tbl_subset(.,"Id",H1) %>% 
  rowwise(Id)%>% mutate(m = mean(c_across(where(is.numeric))))%>%ungroup()%>%filter(m>0)
cage_enh_HMEC_tbl<-cage_enh_tbl %>% cage_tbl_subset(.,"Id",HMEC) %>% 
  rowwise(Id)%>% mutate(m = mean(c_across(where(is.numeric))))%>%ungroup()%>%filter(m>0)
cage_enh_GM12878_tbl<-cage_enh_tbl %>% cage_tbl_subset(.,"Id",GM12878) %>% 
  rowwise(Id)%>% mutate(m = mean(c_across(where(is.numeric))))%>%ungroup()%>%filter(m>0)

cage_enh_H1_tbl<-cage_enh_H1_tbl %>% 
  cage_tbl_coord_build_fn(.,'Id')
cage_enh_HMEC_tbl<-cage_enh_HMEC_tbl %>% 
  cage_tbl_coord_build_fn(.,'Id')
cage_enh_GM12878_tbl<-cage_enh_GM12878_tbl %>% 
  cage_tbl_coord_build_fn(.,'Id')

save(cage_enh_HMEC_tbl,file="~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_enh_tbl.Rda")
save(cage_enh_H1_tbl,file="~/Documents/multires_bhicect/data/epi_data/H1/CAGE/CAGE_enh_tbl.Rda")
save(cage_enh_GM12878_tbl,file="~/Documents/multires_bhicect/data/epi_data/GM12878/CAGE/CAGE_enh_tbl.Rda")

# build union of TSS and Enh CAGE peak object
HMEC_enh_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_enh_tbl.Rda"
HMEC_tss_file<-"./data/CAGE_tss_coord_HMEC_tbl.Rda"

H1_enh_file<-"~/Documents/multires_bhicect/data/epi_data/H1/CAGE/CAGE_enh_tbl.Rda"
H1_tss_file<-"./data/CAGE_tss_coord_H1_tbl.Rda"

GM12878_enh_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CAGE/CAGE_enh_tbl.Rda"
GM12878_tss_file<-"./data/CAGE_tss_coord_GM12878_tbl.Rda"

cage_peak_HMEC_tbl<-tss_enh_tbl_build_fn(HMEC_enh_file,HMEC_tss_file)

cage_peak_H1_tbl<-tss_enh_tbl_build_fn(H1_enh_file,H1_tss_file)

cage_peak_GM12878_tbl<-tss_enh_tbl_build_fn(GM12878_enh_file,GM12878_tss_file)

save(cage_peak_HMEC_tbl,file="./data/CAGE_union_coord_HMEC_tbl")
save(cage_peak_H1_tbl,file="./data/CAGE_union_coord_H1_tbl")
save(cage_peak_GM12878_tbl,file="./data/CAGE_union_coord_GM12878_tbl")
