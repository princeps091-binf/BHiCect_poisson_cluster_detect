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
