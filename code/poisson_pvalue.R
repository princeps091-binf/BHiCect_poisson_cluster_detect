library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
feature_file<-"./data/GRanges/CAGE_union_GM12878_Grange.Rda"
cl_folder<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/BHiCect_Grange/GM12878/"
hic_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
#-----------------------------------------
#Utils. Fn
obj_in_fn<-function(file){
  obj_out<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(obj_out)
}
get_bin_dat<-function(dat_folder,tmp_res,chromo){
  hic_dat<-read_tsv(paste0(hic_dat_folder,tmp_res,"/",chromo,".txt"),col_names = F)
  bin_dat<-tibble(chr=chromo,res=tmp_res,bin=unique(c(hic_dat$X1,hic_dat$X2)))
  return(bin_dat)
}
#-----------------------------------------
feature_GRange<-obj_in_fn(feature_file)

dat_res_set<-grep("b$",list.files(hic_dat_folder),value=T)
cage_chr_reff_l<-do.call(bind_rows,map(dat_res_set,function(tmp_res){
  message(tmp_res)
  tmp_chr_set<-str_split_fixed(grep("^chr",list.files(paste0(hic_dat_folder,tmp_res,"/")),value=T),pattern = "\\.",n=2)[,1]
  tmp_res_tbl<-do.call(bind_rows,lapply(tmp_chr_set,function(chromo){
    bin_dat<-get_bin_dat(hic_dat_folder,tmp_res,chromo)
    bin_Grange<-GRanges(seqnames=bin_dat$chr,
                        ranges = IRanges(start=as.numeric(bin_dat$bin),
                                         end=as.numeric(bin_dat$bin) + res_num[bin_dat$res] -1
                        ))
    return(tibble(chr=chromo,res=tmp_res,chr.cage.count=sum(countOverlaps(bin_Grange,feature_GRange)),bin.count=length(bin_Grange)))
    
  }))
}))

chr_set<-str_split_fixed(grep("^chr",list.files(cl_folder),value=T),"_",n=2)[,1]
cl_res_l<-vector('list',length(chr_set))
names(cl_res_l)<-chr_set
for(chromo in chr_set){
  message(chromo)
  cl_chr_tbl<-obj_in_fn(paste0(cl_folder,chromo,"_BHiCect_cl.Rda"))
  chr_feature_Grange<-feature_GRange[seqnames(feature_GRange)==chromo]
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl,c("chr_feature_Grange"))#
  cl_inter_vec<-unlist(parLapply(cl,cl_chr_tbl$GRange,function(x){
    sum(countOverlaps(x,chr_feature_Grange))
  }))
  stopCluster(cl)
  rm(cl)
  
  cl_chr_tbl<-cl_chr_tbl%>%mutate(feature_n=cl_inter_vec)%>% filter(feature_n>0)
  
  cl_res_l[[chromo]]<-cl_chr_tbl %>% 
    dplyr::select(chr,cl,res,feature_n) %>% 
    mutate(cl.bin=as.numeric(str_split_fixed(cl,"_",3)[,2]))
}
cl_res_tbl<-do.call(bind_rows,cl_res_l)
plan(multisession,workers=5)
cl_res_tbl<-cl_res_tbl %>% 
  left_join(.,cage_chr_reff_l) %>% 
#  dplyr::slice(1:5) %>% 
  mutate(pois.pval=future_pmap_dbl(list(feature_n,cl.bin,chr.cage.count,bin.count),function(feature_n,cl.bin,chr.cage.count,bin.count){
    poisson.test(c(feature_n,chr.cage.count),T=c(cl.bin,bin.count),alternative = "greater")$p.value
  }))
plan(sequential)

cl_res_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(feature_n,-log10(pois.pval)))+
  geom_point()+
  facet_wrap(res~.,scales='free')
save(cl_res_tbl,file="./data/GM12878_pois_pval_tbl.Rda")
#------------------------------------------------------------------
load("~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda")


cl_chr_emp_pval_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(feature_n,-log10(emp.pval)))+
  geom_point()+
  facet_wrap(res~.,scales='free')

cl_chr_emp_pval_tbl %>% 
  dplyr::select(chr,res,cl,emp.pval) %>% 
  left_join(.,cl_res_tbl %>% 
              dplyr::select(chr,res,cl,pois.pval)) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(-log10(emp.pval),-log10(pois.pval)))+
  geom_smooth()+
  facet_wrap(res~.,scales='free')

