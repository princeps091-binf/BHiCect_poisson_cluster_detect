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
cl_folder<-"./data/GRanges/BHiCect_Grange/GM12878/"
hic_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/100kb/"
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
  hic_dat<-obj_in_fn(paste0(hic_dat_folder,tmp_res,"/",chromo,".txt"))
  bin_dat<-tibble(chr=chromo,res=tmp_res,bin=unique(c(hic_dat$X1,hic_dat$X2)))
  return(bin_dat)
}
#-----------------------------------------
